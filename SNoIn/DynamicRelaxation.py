"""
The ``DynamicRelaxation`` module:
-------------------
Using the concept of Dynamic Relaxation, this module provides a solution
for nonlinear problems with large deformation.
"""

from dolfin import *
from tIGAr.timeIntegration import *

def DynamicRelaxationSolve(spline, R, J, u, DENS, damp, max_count=3, tol=1e-3):
    """
    Solves the nonlinear problem using Dynamic Relaxation. "timeInt" is the
    "BackwardEulerIntegrator" class from "tIGAr.timeIntegration" module.
    Inertial and damping contribution will be added to R and J
    "max_count" and "tol" define the convergence criterion.
    """

    # Quantities from the previous time step (Initialized to be used in time integration)
    u_old = Function(spline.V)
    udot_old = Function(spline.V)

    # Create a time integrator for the displacement.
    DELTA_T = Constant(1)
    timeInt = BackwardEulerIntegrator(DELTA_T, u,
                                      (u_old, udot_old))
    # Update DELTA_T
    def updateDt(ratio):
        timeInt.DELTA_T.assign(timeInt.DELTA_T * ratio)
        timeInt.DELTA_T_reciprocal.assign(1 / (timeInt.DELTA_T * ratio))

    # Inertial contribution to the residual:
    z_hom = TestFunction(spline.V)
    acceleration = timeInt.xddot()
    velocity = timeInt.xdot()
    dWmass = DENS * inner(acceleration, z_hom) * spline.dx
    dWdamp = damp * inner(velocity, z_hom) * spline.dx

    timeInt.x_old.assign(u)
    InitialV = Constant(0)  # Initial velocity of the beam
    timeInt.xdot_old.assign(spline.project(Constant([0, 0, InitialV]), rationalize=False))

    R += dWmass + dWdamp
    J += derivative(dWmass + dWdamp, u)

    # Solve
    i = 0
    converge = False
    count = 0  # For convergence criterion
    while not converge:
        print("------- Time step "+str(i+1)
                  +" , t = "+str(timeInt.t)+" -------")
        print('dt =', float(timeInt.DELTA_T))

        u_temp = Function(spline.V)
        try:
            u_temp.assign(u)
            spline.solveNonlinearVariationalProblem(R, J, u)  # Newton iteration

            # Update the quantities of time integration, and double the Delta_T
            timeInt.advance()
            updateDt(2)
            i += 1

        except:
            # If the nonlinear solver doesn't converge, enter here and half the Delta_T and rerun.
            timeInt.updateDt(0.5)
            # For convergence criterion, see if the shell moves less than the given tolerance
            if abs(norm(u)/norm(u_temp)-1) <= tol:
                count += 1
            u.assign(u_temp)

        # Convergence criterion, if beam has not moved for 'count' times, we've converged.
        if count == max_count:
            converge = True
            print('Converged!')
