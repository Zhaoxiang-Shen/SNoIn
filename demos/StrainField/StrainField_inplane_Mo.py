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

# Strain tensor
A0, A1, A2, _, A, B = surfaceGeometry(spline, X)
a0, a1, a2, _, a, b = surfaceGeometry(spline, x)
epsilon = 0.5 * (a - A)
kappa = B - b
epsilonBar = covariantRank2TensorToCartesian2D(epsilon, A, A0, A1)
kappaBar = covariantRank2TensorToCartesian2D(kappa,A,A0,A1)

def eig(M):
    """
    A Closed-form of the eigen decomposition of a 2x2 Hermitian matrix
    https://hal.archives-ouvertes.fr/hal-01501221/document
    """
    delta = sqrt(4*M[1,0]**2+(M[0,0]-M[1,1])**2)
    lambda1 = (M[0,0]+M[1,1]+delta)/2
    lambda2 = (M[0,0]+M[1,1]-delta)/2
    s1 = as_vector([lambda1-M[1,1],M[1,0]])
    s2 = as_vector([lambda2-M[1,1],M[1,0]])

    # In case M[1,0] == 0:
    # s1 = conditional(gt(M[1,0]**2,Constant(1e-6)), unit(s1), as_vector([1,0]))
    # s2 = conditional(gt(M[1,0]**2,Constant(1e-6)), unit(s2), as_vector([0,1]))
    return lambda1,lambda2,unit(s1),unit(s2)

C = Constant(0.1)
k1,k2,N1,N2 = eig(kappaBar)
e_Mo = epsilonBar + C * (k1**2*outer(N1,N1)+k2**2*outer(N2,N2))
e_Mo_max = 0.5 * (e_Mo[0,0]+e_Mo[1,1]) \
             + 0.5 * sqrt((e_Mo[0,0]-e_Mo[1,1])**2 + 4 * e_Mo[0,1]**2)

e_Mo_max = spline_sub.project(e_Mo_max, rationalize=False, lumpMass=True)
e_Mo_11 = spline_sub.project(e_Mo[0,0], rationalize=False, lumpMass=True)
e_Mo_22 = spline_sub.project(e_Mo[1,1], rationalize=False, lumpMass=True)
e_Mo_12 = spline_sub.project(e_Mo[0,1], rationalize=False, lumpMass=True)

###
e_Mo_max.rename("e_Mo_max", "e_Mo_max")
File(folder_name + "/e_Mo_max.pvd") << e_Mo_max

e_Mo_11.rename("e_Mo_11", "e_Mo_11")
File(folder_name + "/e_Mo_11.pvd") << e_Mo_11

e_Mo_22.rename("e_Mo_22", "e_Mo_22")
File(folder_name + "/e_Mo_22.pvd") << e_Mo_22

e_Mo_12.rename("e_Mo_12", "e_Mo_12")
File(folder_name + "/e_Mo_12.pvd") << e_Mo_12
