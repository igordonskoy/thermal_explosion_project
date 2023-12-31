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


def lin_diff(x,f):
    # 1st order approximation of derivative
    N = len(x);
    df = np.zeros((N,1));
    for i in range(0,N-1):
        df[i] = (f[i+1]-f[i])/(x[i+1]-x[i]);
    df[N-1] = df[N-2];
    return df


def quadro_diff(x,f):
# 2nd order Lagrange polynom approximation of derivatives
    N = len(x);
    if len(x) != len(f):
        print('Dimensions should be equal!')
        return 0
    elif len(x) == 1:
        print('Not enough points!')
        return 0
    else:
        df = np.zeros((N,1));
        df[0] = (f[1]-f[0])/(x[1]-x[0]);                  # beginning and ending
        df[N-1] = (f[N-1]-f[N-2])/(x[N-1]-x[N-2]);        # are approximated by straight differences
        for i in range(1,N-1):
            a1 = f[i-1]/(x[i-1]-x[i])/(x[i-1]-x[i+1]);
            a2 = f[i]/(x[i]-x[i-1])/(x[i]-x[i+1]);
            a3 = f[i+1]/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
            df[i] = x[i-1]*(-a2-a3) + x[i+1]*(-a1-a2) + x[i]*(a1 + a3+2*a2);
        return df
