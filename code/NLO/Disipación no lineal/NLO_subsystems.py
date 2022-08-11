from CosimUnit import CosimUnit
import numpy as np

# ________________________________________________________________________________________________YF_US
# This file defines a class that represents a mass-spring system that
# returns force as an output.

#   |         |----|             
#   |--/\/\/--| m1 |----/\/\/---- -> f
#   | k1, c1  |----|   k_c, c_c  

# Input: 	displacement and velocity of the second subsystem
# Output: 	coupling force exerted on environment
# The coupling damping and stiffness are required by this subsystem
class mass_yF_uS(CosimUnit):

    # _______________________________________________________ CONSTRUCTOR
    def __init__(self, name, SYS, idx):
        # Call base class
        super().__init__(name)

        if idx != 1 and idx != 2:
            raise ValueError('Incorrect subsystem index')

        # Attributes
        self.__stepSize = 0
        self.__idx      = idx
        self.__m        = SYS.m[self.__idx - 1]
        self.__k        = SYS.k[self.__idx - 1]
        self.__c        = SYS.c[self.__idx - 1]
        self.__s0       = SYS.s0[self.__idx - 1]
        self.__sd0      = SYS.sd0[self.__idx - 1]
        self.__kc       = SYS.kc
        self.__cc       = SYS.cc

        # Time
        self.__time     = 0.0   # Current time

        # Inputs
        self.__u        = np.zeros(2)

        # States
        self.__x        = self.__s0
        self.__xd       = self.__sd0
        self.__xdd      = 0.0

        # Forces
        self.__fc_sp    = 0.0
        self.__fc_dp    = 0.0
        self.__f_sp     = 0.0
        self.__f_dp     = 0.0

    def evaluateSpringForces(self):

        # Las fuerzas del amortiguador ahora son calculadas como f = c * sign(xd) * xd^2
        # Coupling spring and damper
        self.__fc_sp  = -(self.__kc * (self.__x - self.__u[0]))
        self.__fc_dp  = -(self.__cc * np.sign(self.__xd - self.__u[1]) * ((self.__xd - self.__u[1])**2))

        # Spring and damper betweeen mass and wall
        self.__f_sp = - self.__k * self.__x
        self.__f_dp = - self.__c * np.sign(self.__xd) * (self.__xd**2)
    
    # ___________________________________________________ METHODS, INHERITED
    #
    # These methods are user-defined
    # Co-simulation procedures
    def doStep(self, finalT):
        # Call base class function implementation
        super().doStep(finalT)

        # Additional implementation
        # Take steps until t goes beyond finalT
        while finalT - self.__time > self.__stepSize / 2.0:
            # Evaluate accelerations
            self.evaluateSpringForces()
            LHS = self.__m
            RHS = self.__fc_sp + self.__fc_dp + self.__f_sp + self.__f_dp
            self.__xdd = RHS / LHS

            # Integrate (semi-implicit forward Euler)
            self.__xd    += self.__stepSize * self.__xdd
            self.__x     += self.__stepSize * self.__xd

            # Increase time
            self.__time += self.__stepSize
                
            # Evaluate outputs
            self.evaluateSpringForces()

    def initialize(self, t):
        # Call base class function implementation
        super().initialize(t)

        # Additional implementation
        print('Initializing mass block', self.__idx)

        # Evaluate initial accelerations
        self.evaluateSpringForces()
        LHS = self.__m
        RHS = self.__fc_sp + self.__fc_dp + self.__f_sp + self.__f_dp
        self.__xdd = RHS / LHS

    def terminate(self):
        # Call base class function implementation
        super().terminate()

        # Additional implementation
        print('Terminating mass block ', self.__idx, 'time: ', round(self.__time, 4), 's')

    def setInputs(self, u):
        # Call base class function implementation
        super().setInputs(u)

        # Additional implementation
        # Position and velocity of external subsystem
        self.__u = u
    
    def  readOutputs(self):
        # Call base class function implementation
        super().readOutputs()

        # Additional implementation
        # The output of this subsystem is the coupling force 
        # that the subsystem exerts on its environment
        # We also return position and velocity to have access to all
        # the system information
        return np.array((-(self.__fc_sp + self.__fc_dp), self.__x, self.__xd, self.__xdd))

    # _______________________________________________________ ACCESSORS
    def setStepSize(self, dt):
        self.__stepSize = dt    

# ________________________________________________________________________________________________YS_UF
# This file defines a class that represent a mass-spring system that
# returns displacement and velocity as an output.

#   |         |----|             
#   |--/\/\/--| m1 |----o
#   | k1, c1  |----| 

# Input: coupling force that the subsystem exerts on its environment
# Output: displacement and velocity of the second subsystem
# The coupling damping and stiffness are not included in this subsystem
class mass_yS_uF(CosimUnit):

    # _______________________________________________________ CONSTRUCTOR
    def __init__(self, name, SYS, idx):
        # Call base class
        super().__init__(name)
        
        if idx != 1 and idx != 2:
            raise ValueError('Incorrect subsystem index')

        # Attributes
        self.__stepSize = 0
        self.__idx      = idx
        self.__m        = SYS.m[self.__idx - 1]
        self.__k        = SYS.k[self.__idx - 1]
        self.__c        = SYS.c[self.__idx - 1]
        self.__s0       = SYS.s0[self.__idx - 1]
        self.__sd0      = SYS.sd0[self.__idx - 1]

        # Time
        self.__time     = 0.0   # Current time

        # Inputs
        self.__u        = np.zeros(1)

        # States
        self.__x        = self.__s0
        self.__xd       = self.__sd0
        self.__xdd      = 0.0

        # Forces
        self.__f_sp     = 0.0
        self.__f_dp     = 0.0

    def evaluateSpringForces(self):
        # Las fuerzas del amortiguador ahora son calculadas como f = c * sign(xd) * xd^2
        # Spring and damper betweeen mass and wall
        self.__f_sp = - self.__k * self.__x
        self.__f_dp = - self.__c * np.sign(self.__xd) * (self.__xd**2)

    # ___________________________________________________ METHODS, INHERITED
    def initialize(self, t):
        # Call base class function implementation
        super().initialize(t)

        # Additional implementation
        print('Initializing mass block', self.__idx)

        # Evaluate initial accelerations
        self.evaluateSpringForces()
        LHS = self.__m
        RHS = self.__u[0] + self.__f_sp + self.__f_dp
        self.__xdd = RHS / LHS

    def doStep(self, finalT):
        # Call base class function implementation
        super().doStep(finalT)
        
        # Additional implementation
        # Take steps until t goes beyond finalT
        while finalT - self.__time > self.__stepSize / 2.0:
            # Evaluate accelerations
            self.evaluateSpringForces()
            LHS = self.__m
            RHS = self.__u[0] + self.__f_sp + self.__f_dp
            self.__xdd = RHS / LHS

            # Integrate (semi-implicit forward Euler)
            self.__xd    += self.__stepSize * self.__xdd
            self.__x     += self.__stepSize * self.__xd

            # Increase time
            self.__time += self.__stepSize

    def terminate(self):
        # Call base class function implementation
        super().terminate()

        # Additional implementation
        print('Terminating mass block ', self.__idx, 'time: ', round(self.__time, 4), 's')

    def setInputs(self, u):
        # Call base class function implementation
        super().setInputs(u)

        # Additional implementation
        # Coupling force 
        # Need to change sign: from "exerted on the environment" to
        # "exerted on the system"
        self.__u = -u

    def  readOutputs(self):
        # Call base class function implementation
        super().readOutputs()

        # Additional implementation
        # The output of this subsystem is the position and velocity of the mass
        return np.array((self.__x, self.__xd, self.__xdd))

    # _______________________________________________________ ACCESSORS
    def setStepSize(self, dt):
        self.__stepSize = dt  
