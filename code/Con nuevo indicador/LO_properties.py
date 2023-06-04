import numpy as np

# Oscillator properties
class sys:

    def __init__(self, m1, m2, k1, k2, kc, c1, c2, cc, s10, s1d0, s20, s2d0):
        # System matrices
        self.M          = np.zeros((2, 2))
        self.M[0, 0]    = m1
        self.M[1, 1]    = m2

        self.K          = np.zeros((2, 2))
        self.K[0, 0]    = k1 + kc
        self.K[0, 1]    = - kc
        self.K[1, 0]    = - kc
        self.K[1, 1]    = k2 + kc

        self.C          = np.zeros((2, 2))
        self.C[0, 0]    = c1 + cc
        self.C[0, 1]    = - cc
        self.C[1, 0]    = - cc
        self.C[1, 1]    = c2 + cc

        # System parameters
        self.m          = np.array((m1, m2))
        self.k          = np.array((k1, k2))
        self.c          = np.array((c1, c2))
        self.kc         = kc
        self.cc         = cc
        self.s0         = np.array((s10, s20))
        self.sd0        = np.array((s1d0, s2d0))