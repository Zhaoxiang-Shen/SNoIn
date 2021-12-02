"""
The ``LJpotential`` module:
-------------------
Module for LJ potential (vdW) between a MoS2 molecule and a Si3N4 half space
The Si3N4 half space is approximated as homogeneous solid.
"""

from dolfin import *

# Material parameters, rho(atom/nm**3),eps(eV),sigma(nm):
ev2kgm2 = 1.60218e-19
ev2kgnm2 = ev2kgm2*1e18
rho_Si = Constant(41.8552)
eps_Mo_Si = Constant(0.038120 * ev2kgnm2)
sigma_Mo_Si = Constant(0.33022)
eps_S_Si = Constant(0.001799*ev2kgnm2)
sigma_S_Si = Constant(0.37106)

rho_N = Constant(55.8070)
eps_Mo_N = Constant(0.077892*ev2kgnm2)
sigma_Mo_N = Constant(0.30261)
eps_S_N = Constant(0.003672*ev2kgnm2)
sigma_S_N = Constant(0.34345)

# Atom densities for a MoS2 layer
rho_Mo = Constant(11.35)   # (atom/nm**2)
rho_S = Constant(11.35)   # (atom/nm**2)

# Vertical distance between Mo and S atoms.
d_S1_Mo = Constant(0.1595)
d_S2_Mo = Constant(-0.1595)

# Pairwise LJ potential, depending on distance squared
def phi(r2, rho, eps, sigma):
    sr6 = (sigma*sigma/r2)**3
    return 4.0*rho*eps*(sr6**2 - sr6)

# Semi-infinite integration of 'phi()', as the LJ potential between a particle and a half space
def I_3D_sub(h,rho,eps,sigma):
    gamma = 6/5 * pi * eps * rho * sigma ** 2
    h0 = sigma
    return - gamma * (5/3 * h0**4/(3*h**3) - 2/3 * h0**10/(9*h**9))

# LJ energy density of a monolayer MoS2.
# Summing the contribution from all different atom pairs, i.e. Mo-Si/N, S-Si/N.
def I_MoS2_Si3N4(d_Mo, n_a2=1):
    """
    For a MoS2 molecule, “d_S == d_Mo + d_S_Mo * n_a2”,
    where n_a2 == abs(unit(inner(n,a2))),
    “n” is the normal vector of the tangent plane, "a2” is the normal vector at "x"
    (In the demos, "n_a2=1" can produce results with no apparent difference).
    """
    shift_S1_Mo = d_S1_Mo * n_a2
    shift_S2_Mo = d_S2_Mo * n_a2
    d_S1 = d_Mo + shift_S1_Mo
    d_S2 = d_Mo + shift_S2_Mo
    I_sum = rho_Mo * (I_3D_sub(d_Mo,rho_N,eps_Mo_N,sigma_Mo_N) + I_3D_sub(d_Mo,rho_Si,eps_Mo_Si,sigma_Mo_Si)) + \
            rho_S * (I_3D_sub(d_S1,rho_N,eps_S_N,sigma_S_N) + I_3D_sub(d_S1,rho_Si,eps_S_Si,sigma_S_Si) + \
                      I_3D_sub(d_S2,rho_N,eps_S_N,sigma_S_N) + I_3D_sub(d_S2,rho_Si,eps_S_Si,sigma_S_Si))
    return I_sum
