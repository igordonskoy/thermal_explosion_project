import numpy as np


def transfer_matrix(N, lbc=1):
    A = np.zeros((N, N))
    for i in range(1,N-1):
        A[i,i] = -2
        A[i,i-1] = 1
        A[i,i+1] = 1
    A[0,0] = 1
    if lbc == 2:
        A[0,1] = -1
    A[N-1,N-1] = 1
    return A


def progonka(A, B): # TDMA realization 
    # A and B are numpy arrays of NxN abd Nx1 shape correspondingly
    c = np.zeros((len(B), 1))  # empty vessels
    d = np.zeros((len(B), 1))  # for coefficients
    X = np.zeros((len(B), 1))  # and for solution
    c[0] = A[0,1] / A[0,0]
    d[0] = B[0] / A[0,0]
    for i in range(1, len(B)-1):       # forward cycle
        c[i] = A[i,i+1] / (- A[i,i-1] * c[i-1] + A[i,i])
        d[i] = (B[i] - A[i,i-1]*d[i-1]) / (- A[i,i-1] * c[i-1] + A[i,i])
    X[-1] = (B[-1] - A[-1, -2]*d[-2]) / (A[-1,-1] -A[-1,-2]*c[-2] )
    for i in range(len(B)-2, -1, -1):  # backward cycle
        X[i] = - c[i]*X[i+1] + d[i]
    return X

