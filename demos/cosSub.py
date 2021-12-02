from tIGAr import *
from tIGAr.NURBS import *
from tIGAr.BSplines import *
from numpy import array
from igakit.nurbs import NURBS as NURBS_ik
from igakit.io import PetIGA
from ShNAPr.kinematics import *
from ShNAPr.hyperelastic import *
from ShNAPr.SVK import *
from SNoIn.Substrate import *
from SNoIn.DynamicRelaxation import *
from SNoIn.LJpotential import *
from SNoIn.TangentMethod import *

import time
time_start = time.time()

# Use TSFC representation, due to complicated forms:
parameters["form_compiler"]["representation"] = "tsfc"
import sys
sys.setrecursionlimit(1000000)

####### Spline setup for IGA #######

# Polynomial degree of the basis functions; must be >1, because
# functions need to be at least $C^1$ for the displacement-only 4-th-order
# beam formulation.
p = 2

# Number of elements in each direction:
Nel = 25

# Length of the shell:
scale_geo = 2
L = 20 * scale_geo
Lu = L
Lv = L*sin(pi/3)
print('L =', L, ', Nel=', Nel)

# Open knot vectors for a_dom one-Bezier-element bi-unit square.
uKnots = [0.0, 0.0, 1.0, 1.0]
vKnots = [0.0, 0.0, 1.0, 1.0]

# Array of control points, define the shape.
cpArray = array([[[0.0, 0.0], [Lu, 0.0]],
                 [[-L*sin(pi/6), Lv], [L*sin(pi/6), Lv]]])
center_coord = (cpArray[0, 0]+cpArray[1, 1])/2

# Create initial mesh
ikNURBS = NURBS_ik([uKnots, vKnots], cpArray)
ikNURBS.elevate(0)
ikNURBS.elevate(1)

# uniform Refinement
numNewKnots = Nel
h = 1/float(numNewKnots)
numNewKnots -= 1
knotList = []
for i in range(0,numNewKnots):
    knotList += [float(i + 1) * h, ]
newKnots = array(knotList)
ikNURBS.refine(0, newKnots)
ikNURBS.refine(1, newKnots)

# Output in PetIGA format
if (mpirank == 0):
    PetIGA().write("out.dat", ikNURBS)
MPI.barrier(worldcomm)

# Degree
degs = [p, p]

# No repetition for periodic BC
kvecs_periodic = [list(ikNURBS.knots[0])[p:-p],
                  list(ikNURBS.knots[1])[p:-p]]

controlMesh = NURBSControlMesh("out.dat", useRect=True,overRefine=0)
field = BSpline(degs, kvecs_periodic)
splineGenerator = FieldListSpline(controlMesh, 3*[field,])

####### Analysis #######
QUAD_DEG = 2*p
spline = ExtractedSpline(splineGenerator, QUAD_DEG)

# substrate
shift_sub = Constant(0)
SubstrateSurface = SubstrateLibrary(spline, controlMesh, field, QUAD_DEG)
A_sub = 16
R_sub = 16
alpha = 2
U_sub = SubstrateSurface.cosSub(center_coord,R_sub,A_sub,alpha,shift=shift_sub)
spline_sub = SubstrateSurface.spline

# Unknown midsurface displacement
y_hom = Function(spline.V)  # in homogeneous representation
y = spline.rationalize(y_hom)  # in physical coordinates
y_hom.assign(spline.project(Constant([0, 0, 0.432]), rationalize=False))

# Reference configuration:
X = spline.F

# Current configuration:
x = X + y

# Elastic constants
# REBO: a: 0.4671307/30.0058e1, b: 0.4854564/26.77892e1
# SW: a: 0.4014012/24.13545e1, b: 0.4098362/22.67569e1

# Shell thickness:
Et_choice = [0,0]  # D=(8.75(SW),17.467(REBO))eV
h_th_list = [[0.4014012,0.4098362], [0.4671307,0.4854564]]  # nm
h_th = Constant(h_th_list[Et_choice[0]][Et_choice[1]])

# The Young's modulus and Poisson ratio:
E_list = [[24.13545e1,22.67569e1], [30.0058e1,26.77892e1]]
nu_list = [0.2686,0.2962]
E = Constant(E_list[Et_choice[0]][Et_choice[1]])  # (kg/(nm*s^2))
nu = Constant(nu_list[Et_choice[0]])

# Elastic energy:
energy_E = surfaceEnergyDensitySVK(spline,X,x,E,nu,h_th)*spline.dx

# LJ potential with analytical I (semi-infinite with no hole):
_, _, a2, _, _, _ = surfaceGeometry(spline, x)  # a2 is the normal vector of shell
d_Mo, n_a2 = d_tangent_plane(spline_sub, U_sub, x[2], y[0], y[1], shift_sub, a2)
energy_LJ = I_MoS2_Si3N4(d_Mo, n_a2) * spline.dx

# Total energy:
energy = energy_E + energy_LJ
R = derivative(energy,y_hom)
J = derivative(R,y_hom)

# Solve
DENS = Constant(1000)
damp = Constant(200)
spline.setSolverOptions(maxIters=10,relativeTolerance=1e-5,linearSolver=PETScLUSolver("mumps"))
DynamicRelaxationSolve(spline, R, J, y_hom, DENS, damp, max_count=3, tol=1e-3)

time_end = time.time()
print('total cost', time_end-time_start)

# Create file for result
d0File = File("result/disp-x.pvd")
d1File = File("result/disp-y.pvd")
d2File = File("result/disp-z.pvd")

(d0, d1, d2) = y_hom.split()
d0.rename("d0", "d0")
d1.rename("d1", "d1")
d2.rename("d2", "d2")
d0File << d0
d1File << d1
d2File << d2

U_sub.rename("U_sub", "U_sub")
File("result/U_sub.pvd") << U_sub
# Geometry files
nsd = 3
for i in range(0, nsd + 1):
    name = "F" + str(i)
    spline.cpFuncs[i].rename(name, name)
    File("result/" + name + "-file.pvd") << spline.cpFuncs[i]

ev2kgm2 = 1.60218e-19
ev2kgnm2 = ev2kgm2*1e18
print('L =', float(L), ', Nel =', Nel, ', Elastic Constant =', Et_choice,'D =', float(1e-18*E*h_th**3/(12*(1-nu**2))/ev2kgm2))
print('A =', float(A_sub), ', r_nonzero =', float(R_sub))
