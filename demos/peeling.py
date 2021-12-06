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
# shell formulation.
p = 2

# Number of elements in each direction:
Nel = 100

# Length of the shell:
L = 34.28
W2L = 20
W = (L/W2L)

# Specify degree in each direction.
degs = [p, p]

# Generate open knot vectors for each direction. Only include periodic BCs in the field
kvecs = [uniformKnots(degs[0],0,L,Nel),
         uniformKnots(degs[1],0,W,int(Nel/W2L))]
kvecs_Pbc = [uniformKnots(degs[0],0,L,Nel,periodic=False),
         uniformKnots(degs[1],0,W,int(Nel/W2L),periodic=True)]
controlMesh = ExplicitBSplineControlMesh(degs,kvecs,extraDim=1)

field = BSpline(degs, kvecs_Pbc)
splineGenerator = FieldListSpline(controlMesh, 3*[field,])

# # Apply Dirichlet BCs to the left, d2 is fixed
parametricDirection = 0
for field_dof in [2]:
    side = 0
    scalarSpline = splineGenerator.getScalarSpline(field_dof)
    sideDofs = scalarSpline.getSideDofs(parametricDirection, side, nLayers=1)
    splineGenerator.addZeroDofs(field_dof, sideDofs)

# # Apply Dirichlet BCs to the right, d0 is fixed
parametricDirection = 0
for field_dof in [0]:
    side = 1
    scalarSpline = splineGenerator.getScalarSpline(field_dof)
    sideDofs = scalarSpline.getSideDofs(parametricDirection, side, nLayers=1)
    splineGenerator.addZeroDofs(field_dof, sideDofs)

####### Analysis #######
QUAD_DEG = 2*p
spline = ExtractedSpline(splineGenerator, QUAD_DEG)

# substrate
shift_sub = Constant(0)
SubstrateSurface = SubstrateLibrary(spline, controlMesh, field, QUAD_DEG)
U_sub = SubstrateSurface.plainSub(shift=shift_sub)
spline_sub = SubstrateSurface.spline

# Unknown midsurface displacement
y_hom = Function(spline.V)  # in homogeneous representation
y = spline.rationalize(y_hom)  # in physical coordinates

# Reference configuration:
X = spline.F

# Current configuration:
x = X + y

# Elastic constants
# REBO: a: 0.4671307/30.0058e1, b: 0.4854564/26.77892e1
# SW: a: 0.4014012/24.13545e1, b: 0.4098362/22.67569e1

# Shell thickness:
Et_choice = [1,0]  # D=(8.75(SW),17.467(REBO))eV
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
# Since the substrate is flat, i.e. n_sub = [0,0,1], it's equivalent to set shift_S_Mo = d_S_Mo * a2[2]
_, _, a2, _, _, _ = surfaceGeometry(spline, x)  # a2 is the normal vector of shell
d_Mo, n_a2 = d_tangent_plane(spline_sub, U_sub, x[2], y[0], y[1], shift_sub, a2)
energy_LJ = I_MoS2_Si3N4(d_Mo, n_a2) * spline.dx

# Total energy:
energy = energy_E + energy_LJ
R = derivative(energy,y_hom)
J = derivative(R,y_hom)

# For simplicity, the different peeling heights are enforced by shifting down the flat substrate,
# while keeping the vertical position of the peeled end (left end).
# In case of high peeling BCs, the 'height' increases gradually,
# and each solution is set to be the initial position for the next 'height'
DENS = Constant(1000)
damp = Constant(200)
spline.setSolverOptions(maxIters=15,relativeTolerance=1e-5,linearSolver=PETScLUSolver("mumps"))
for h_step in [0.5,1,3,5,7,9,11]:
    shift_sub.assign(h_step)
    U_sub.assign(spline_sub.project(Constant(-shift_sub), rationalize=False))
    DynamicRelaxationSolve(spline, R, J, y_hom, DENS, damp, max_count=3, tol=1e-3)

    # Create file for result
    d0File = File("result/peeling_%.2fnm/disp-x.pvd" % float(shift_sub))
    d1File = File("result/peeling_%.2fnm/disp-y.pvd" % float(shift_sub))
    d2File = File("result/peeling_%.2fnm/disp-z.pvd" % float(shift_sub))

    (d0, d1, d2) = y_hom.split()
    d0.rename("d0", "d0")
    d1.rename("d1", "d1")
    d2.rename("d2", "d2")
    d0File << d0
    d1File << d1
    d2File << d2

    U_sub.rename("U_sub", "U_sub")
    File("result/peeling_%.2fnm/U_sub.pvd" % float(shift_sub)) << U_sub

    # Geometry files
    nsd = 3
    for i in range(0, nsd + 1):
        name = "F" + str(i)
        spline.cpFuncs[i].rename(name, name)
        File("result/peeling_%.2fnm/" % float(shift_sub) + name + "-file.pvd") << spline.cpFuncs[i]

time_end = time.time()
print('total cost', time_end-time_start)

ev2kgm2 = 1.60218e-19
ev2kgnm2 = ev2kgm2*1e18
print('L =', float(L), ', Nel =', Nel, ', Elastic Constant =', Et_choice,'D =', float(1e-18*E*h_th**3/(12*(1-nu**2))/ev2kgm2))
