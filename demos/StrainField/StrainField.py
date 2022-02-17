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

# spline.project
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
scale_geo = 1
L = 40 * scale_geo * 1
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
A_sub = 16 * scale_geo
R_sub = 16 * scale_geo
alpha = 2
U_sub = SubstrateSurface.cosSub(center_coord,R_sub,A_sub,alpha,shift=shift_sub)
spline_sub = SubstrateSurface.spline

# Elastic constants
# REBO: a: 0.4671307/30.0058e1, b: 0.4854564/26.77892e1
# SW: a: 0.4014012/24.13545e1, b: 0.4098362/22.67569e1
Et_choice = [1,0]  # D=(8.75(SW),17.467(REBO))eV
P_option = ["SW", "REBO"]
P_choice = P_option[Et_choice[0]]
print(P_choice)
folder_name = "result/cos_L%d_A%d_R%d_Nel%d_" % (int(L),int(A_sub),int(R_sub),int(Nel)) + P_choice

# Unknown midsurface displacement
y_hom = Function(spline.V, folder_name + "/y_hom.xml")
y = spline.rationalize(y_hom)  # in physical coordinates

# Reference configuration:
X = spline.F

# Current configuration:
x = X + y

# Principal Strain
A0, A1, A2, _, A, B = surfaceGeometry(spline, X)
a0, a1, a2, _, a, b = surfaceGeometry(spline, x)
epsilon = 0.5 * (a - A)
epsilonBar = covariantRank2TensorToCartesian2D(epsilon, A, A0, A1)
e_max = 0.5 * (epsilonBar[0,0]+epsilonBar[1,1]) \
             + 0.5 * sqrt((epsilonBar[0,0]-epsilonBar[1,1])**2 + 4 * epsilonBar[0,1]**2)

F_DG = FunctionSpace(spline.mesh, 'DG', 0)
e_max_DG = project(e_max, F_DG)
e_max = spline_sub.project(e_max, rationalize=False, lumpMass=True)
e11 = spline_sub.project(epsilonBar[0,0], rationalize=False, lumpMass=True)
e22 = spline_sub.project(epsilonBar[1,1], rationalize=False, lumpMass=True)
e12 = spline_sub.project(epsilonBar[0,1], rationalize=False, lumpMass=True)

###
e_max.rename("e_max", "e_max")
File(folder_name + "/e_max.pvd") << e_max

e_max_DG.rename("e_max_DG", "e_max_DG")
File(folder_name + "/e_max_DG.pvd") << e_max_DG

e11.rename("e11", "e11")
File(folder_name + "/e11.pvd") << e11

e22.rename("e22", "e22")
File(folder_name + "/e22.pvd") << e22

e12.rename("e12", "e12")
File(folder_name + "/e12.pvd") << e12
