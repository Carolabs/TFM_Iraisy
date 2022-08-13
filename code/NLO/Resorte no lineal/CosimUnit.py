from abc import abstractmethod

# Definition of class CosimUnit
# Sub-system standard interface for co-simulation
# This class is abstract. It cannot be instantiated.
# It is meant to be used as base class for those classes that describe
# actual implementations of co-simulation units.

# The class derives from handle (so we can modify its contents) and from
# matlab.mixin.Heterogeneous, so that CosimMaster can have an array of
# CosimUnits that represent actually different models
class CosimUnit:
    
    # ___________________________________________________ CONSTRUCTOR
    def __init__(self, name):
        # Attributes
        self.__name           = name  # Unit 

    # ___________________________________________________ METHODS, ABSTRACT
    #
    # These methods are user-defined, must be implemented in the classes
    # that inherit from CosimUnit.
    # Co-simulation procedures
    @abstractmethod
    def doStep(self, finalT):
        pass 
    
    @abstractmethod
    def initialize(self, t):
        pass

    @abstractmethod
    def terminate(self):
        pass

    @abstractmethod
    def readOutputs(self):
        pass

    @abstractmethod
    def setInputs(self, u):
        pass