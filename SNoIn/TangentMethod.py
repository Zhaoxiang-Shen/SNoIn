"""
The ``TangentMethod`` module:
-------------------
For an arbitrary point on a shell, find the distance to the tangent plane of its projected point
on the substrate which must defined by the same mesh as the membrane.
"""

from dolfin import *

def unit(v):
    """
    Normalize the vector ``v``.
    """
    return v/sqrt(inner(v,v))

def d_tangent_plane(spline_sub, u_sub, y_shell, ux_shell, uz_shell, a2):
    """
    "u_sub" is the displacement field of the substrate surface,
    "y_shell" is in the current configuration, "ux/uz_shell" is the displacement,
    "a2" is the normal vector at "x".
    """
    # vector
    ux = spline_sub.grad(u_sub)[0]
    uz = spline_sub.grad(u_sub)[1]
    n_u_sub = as_vector([ux, uz, Constant(-1)])  # Normal vector of the tangent plane
    y_vec = as_vector([0,0,Constant(1)])  # A vertical vector
    # angle
    cos_theta = inner(n_u_sub, y_vec) / (sqrt(inner(n_u_sub, n_u_sub)) * sqrt(inner(y_vec, y_vec)))
    cos_theta = sqrt(cos_theta**2)
    # shift coefficient for S atoms
    n_a2 = inner(unit(n_u_sub), a2)
    # distance
    d_Mo = (y_shell - (u_sub + ux * ux_shell + uz * uz_shell)) * cos_theta
    return d_Mo, sqrt(n_a2**2)
