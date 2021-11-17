"""
The ``TangentMethod`` module:
-------------------
For an arbitrary point on a shell, find the distance to the tangent plane of its projected point
on the substrate which must defined by the same mesh as the membrane.
"""

from dolfin import *

def d_tangent_plane(spline_sub, u_sub, y_shell, ux_shell, uz_shell, shift_sub):
    """
    "u_sub" is the displacement field of the substrate surface.
    "y_shell" is in the current configuration, "ux/uz_shell" is the displacement.
    "shift_sub" is the downward shift of the substrate.
    """
    # vector
    ux = spline_sub.grad(u_sub)[0]
    uz = spline_sub.grad(u_sub)[1]
    normal_u_sub = as_vector([ux,uz,Constant(-1)])
    y_vec = as_vector([0,0,Constant(1)])  # A vertical vector
    # angle
    cos_theta = inner(normal_u_sub,y_vec)/(sqrt(inner(normal_u_sub,normal_u_sub))*sqrt(inner(y_vec,y_vec)))
    cos_theta = sqrt(cos_theta**2)
    # distance
    d = (y_shell - (u_sub + ux * ux_shell + uz * uz_shell - shift_sub)) * cos_theta
    return d
