"""
The ``timeIntegration`` module
(modified module, the tIGAr version is not compatible for SNoIn)
------------------------------
is not strictly related to IGA, but contains routines for time integration
methods commonly used in conjunction with IGA, that we found convenient for 
implementing demos.
"""

from tIGAr.common import *

class BackwardEulerIntegrator:

    """
    Class to encapsulate backward Euler formulas for first- and second-order
    ODE systems.  
    """

    def __init__(self,DELTA_T,x,oldFunctions,t=0.0):
        """
        Initialize a backward Euler integrator with time step ``DELTA_T``.
        The unknown function is ``x``.  The sequence of ``Function``
        objects ``oldFunctions`` provides data from the previous time step.
        If ``oldFunctions`` contains only one ``Function``, then the system
        is assumed to be of order 1 and that function is interpreted as the
        initial value of ``x``.  If ``oldFunctions`` contains an additional
        element, then the ODE system is assumed to be of second order, and
        this additional element is interpreted as the initial velocity.
        The parameter ``t`` is the initial time, and defaults to zero.
        """
        self.systemOrder = len(oldFunctions)
        self.DELTA_T = DELTA_T
        self.DELTA_T_reciprocal = Constant(1.0/DELTA_T)
        self.x = x
        self.x_old = oldFunctions[0]
        if(self.systemOrder == 2):
            self.xdot_old = oldFunctions[1]
        self.t = t + float(DELTA_T) # DELTA_T may be a Constant already
            
    def xdot(self):
        """
        Returns the approximation of the velocity at the current time step.
        """

        return self.DELTA_T_reciprocal*self.x \
            - self.DELTA_T_reciprocal*self.x_old

    def xddot(self):
        """
        Returns the approximation of the acceleration at the current time
        step.
        """

        return self.DELTA_T_reciprocal * self.xdot() \
           - self.DELTA_T_reciprocal * self.xdot_old

            
    def advance(self):
        """
        Overwrites the data from the previous time step with the
        data from the current time step.
        """
        x_old = Function(self.x.function_space())
        x_old.assign(self.x)
        if(self.systemOrder==2):
            xdot_old = Function(self.x.function_space())
            xdot_old.assign(self.xdot())
        self.x_old.assign(x_old)
        if(self.systemOrder==2):
            self.xdot_old.assign(xdot_old)
        self.t += float(self.DELTA_T)

    def updateDt(self,ratio):
        self.DELTA_T.assign(self.DELTA_T * ratio)
        self.DELTA_T_reciprocal.assign(1 / (self.DELTA_T * ratio))
