import numpy as np
import scipy as sp
from scipy.linalg import expm

# Function to solve the LO
def analyticalSolution(sys, H, t):

    # Parameters
    m1 = sys.m[0]
    m2 = sys.m[1]
    k1 = sys.k[0]
    k2 = sys.k[1]
    kc = sys.kc
    c1 = sys.c[0]
    c2 = sys.c[1]
    cc = sys.cc
    s10 = sys.s0[0]
    s20 = sys.s0[1]
    s1d0 = sys.sd0[0]
    s2d0 = sys.sd0[1]

    # Number of timesteps
    nsteps = round(t / H) + 1

    # Lead matrix of the system
    L       = np.zeros((4, 4))
    L[0, 1] = 1.0
    L[1, 0] = -(k1 + kc) / m1
    L[1, 1] = -(c1 + cc) / m1
    L[1, 2] = kc / m1
    L[1, 3] = cc / m1

    L[2, 3] = 1.0
    L[3, 0] = kc / m2
    L[3, 1] = cc / m2
    L[3, 2] = -(k2 + kc) / m2
    L[3, 3] = -(c2 + cc) / m2

    # Initial conditions vector
    s0    = np.zeros((4, ))
    s0[0] = s10
    s0[1] = s1d0
    s0[2] = s20
    s0[3] = s2d0

    # Mass and stiffness matrices (these do not include coupling)
    M       = np.zeros((2, 2))
    M[0, 0] = m1
    M[1, 1] = m2

    K       = np.zeros((2, 2))
    K[0, 0] = k1 + kc
    K[0, 1] = - kc
    K[1, 0] = - kc
    K[1, 1] = k2 + kc

    K1       = np.zeros((2, 2))
    K1[0, 0] = k1
    K1[1, 1] = k2

    C       = np.zeros((2, 2))
    C[0, 0] = c1 + cc
    C[0, 1] = - cc
    C[1, 0] = - cc
    C[1, 1] = c2 + cc

    C1       = np.zeros((2, 2))
    C1[0, 0] = c1
    C1[1, 1] = c2
    
    # Eigenvalues and vectors
    A = np.vstack((np.hstack((M, np.zeros((2, 2)))), np.hstack((np.zeros((2, 2)), np.identity(2)))))
    B = np.vstack((np.hstack((C, K)), np.hstack((-np.identity(2), np.zeros((2, 2))))))
    eigenValues, eigenVectors = np.linalg.eig(-np.linalg.solve(A, B))

    # Numeric evaluation of the solution
    t        = np.zeros(nsteps)
    s1       = np.zeros(nsteps)
    s1d      = np.zeros(nsteps)
    s1dd     = np.zeros(nsteps)
    s2       = np.zeros(nsteps)
    s2d      = np.zeros(nsteps)
    s2dd     = np.zeros(nsteps)
    T        = np.zeros(nsteps)
    V        = np.zeros(nsteps)
    E        = np.zeros(nsteps)
    Wnc      = np.zeros(nsteps)
    WncAccSt = np.zeros(nsteps)
    Fc       = np.zeros(nsteps)

    # Current step
    step        = 0
    s1[0]       = s10
    s1d[0]      = s1d0
    s1dd[0]     = 0
    s2[0]       = s20
    s2d[0]      = s2d0
    s2dd[0]     = 0

    # Accumulated work of non-conservative forces
    WncAcc = 0.0

    # Compute system dynamics
    while step < nsteps:
        # Update time
        if step > 0:
            t[step] = t[step - 1] + H

        # Evaluate kinematics
        exp = expm(L * t[step])
        z = np.matmul(exp, s0)
        s1[step]    = z[0]
        s1d[step]   = z[1]
        s2[step]    = z[2]
        s2d[step]   = z[3]

        # Evaluate dynamics
        Q           = - np.matmul(K, np.array((s1[step], s2[step]))) - np.matmul(C, np.array((s1d[step], s2d[step])))
        acc         = np.linalg.solve(M, Q)
        s1dd[step]  = acc[0]
        s2dd[step]  = acc[1]

        # Evaluate coupling
        Fc[step]    = kc * (s2[step] - s1[step]) + cc * (s2d[step] - s1d[step])
        
        # Evaluate energy
        T[step]     = 0.5 * np.matmul(np.matmul(np.array((s1d[step], s2d[step])).T, M), np.array((s1d[step], s2d[step])))
        V[step]     = 0.5 * np.matmul(np.matmul(np.array((s1[step], s2[step])).T, K), np.array((s1[step], s2[step])))
        E[step]     = T[step] + V[step]

        # Nonconservative work
        Wnc[step]       = np.matmul(-(-np.matmul(C1, np.array((s1d[step], s2d[step]))) + cc * (s2d[step] - s1d[step]) * np.array((1, -1))), np.array((s1d[step], s2d[step]))) * H
        WncAcc          += Wnc[step]
        WncAccSt[step]  = WncAcc

        # Update step
        step += 1

    return t, s1, s1d, s1dd, s2, s2d, s2dd, Fc, T, V, E, Wnc, WncAccSt, eigenValues, eigenVectors