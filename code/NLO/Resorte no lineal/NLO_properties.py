import numpy as np

# Oscillator properties
class sys:

    def __init__(self, m1, m2, k1, k2, kc, c1, c2, cc, s10, s1d0, s20, s2d0):

        # System parameters
        self.m          = np.array((m1, m2))
        self.k          = np.array((k1, k2))
        self.c          = np.array((c1, c2))
        self.kc         = kc
        self.cc         = cc
        self.s0         = np.array((s10, s20))
        self.sd0        = np.array((s1d0, s2d0))