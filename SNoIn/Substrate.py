from dolfin import *
from tIGAr import *
from tIGAr.NURBS import *
from tIGAr.BSplines import *

def dist(x, z, center):
    d = ((x - center[0]) ** 2 + (z - center[1]) ** 2) ** 0.5
    return d

class SubstrateLibrary:

    def __init__(self,spline_shell,ControlMesh,field,deg):
        self.splineGenerator = FieldListSpline(ControlMesh, [field,])

        # Write extraction data to the filesystem.
        DIR = "./extraction"
        self.splineGenerator.writeExtraction(DIR)
        self.spline = ExtractedSpline(DIR, deg, mesh=spline_shell.mesh)
        self.x = self.spline.spatialCoordinates()

    def plainSub(self, shift=0):
        u_sub = self.spline.project(Constant(-shift), rationalize=False)
        return u_sub

    def trenchSub(self,center,R,A,alpha,shift=0):
        zero_region = conditional(gt(((self.x[0]-center[0])**2)**0.5, Constant(R)), 0, Constant(1))
        wr = Constant(1/R)
        cos_r = (cos(pi * wr * (self.x[0] - center[0])) + 1) / 2
        u_sub = self.spline.project(zero_region * Constant(-A) * cos_r ** alpha
                                        + Constant(-shift), rationalize=False)
        return u_sub

    def cosSub(self,center,R,A,alpha,shift=0):
        d = dist(self.x[0], self.x[1], center)
        zero_region = conditional(gt(d, Constant(R)), 0, Constant(1))
        wr = Constant(1/R)
        cos_r = (cos(pi * wr * ((self.x[0] - center[0]) ** 2 + (self.x[1] - center[1]) ** 2) ** 0.5) + 1)/2
        u_sub = self.spline.project(zero_region * Constant(-A) * cos_r ** alpha
                                        + Constant(-shift), rationalize=False)
        return u_sub

