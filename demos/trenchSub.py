from tIGAr import *
from tIGAr.NURBS import *
from tIGAr.BSplines import *
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
Nel = 50

# Length of the shell:
scale_geo = 2
L = 20 * scale_geo
W2L = 15
W = int(L/W2L)
print('L =', L, ', Nel=', Nel)
center_coord = [L/2, W/2]

# Specify degree in each direction.
degs = [p, p]

# Generate open knot vectors for each direction. Only include periodic BCs in the field
kvecs = [uniformKnots(degs[0],0,L,Nel),
         uniformKnots(degs[1],0,W,int(Nel/W2L))]
kvecs_Pbc = [uniformKnots(degs[0],0,L,Nel,periodic=True),
         uniformKnots(degs[1],0,W,int(Nel/W2L),periodic=True)]
controlMesh = ExplicitBSplineControlMesh(degs,kvecs,extraDim=1)

field = BSpline(degs, kvecs_Pbc)
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
U_sub = SubstrateSurface.trenchSub(center_coord,R_sub,A_sub,alpha,shift=shift_sub)
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
Constant_choice = 0  # D=(8.75(SW),17.467(REBO))eV
h_th_list = [0.4014012, 0.4671307]  # nm
h_th = Constant(h_th_list[Constant_choice])

# The Young's modulus and Poisson ratio:
E_list = [24.13545e1, 30.0058e1]
nu_list = [0.2686,0.2962]
E = Constant(E_list[Constant_choice])   # (kg/(nm*s^2))
nu = Constant(nu_list[Constant_choice])

# Elastic energy:
energy_E = surfaceEnergyDensitySVK(spline,X,x,E,nu,h_th)*spline.dx

# LJ potential with analytical I (semi-infinite with no hole):
d_Mo = d_tangent_plane(spline_sub, U_sub, x[2], y[0], y[1], shift_sub)
energy_LJ = I_MoS2_Si3N4(d_Mo) * spline.dx

# Energy summation
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
print('L =', float(L), ', Nel =', Nel, ', Elastic Constant =', Constant_choice,'D =', float(1e-18*E*h_th**3/(12*(1-nu**2))/ev2kgm2))
print('A =', float(A_sub), ', r_nonzero =', float(R_sub))
