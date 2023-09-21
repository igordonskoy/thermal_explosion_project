def progonka_col(A, B):
    c = np.zeros((len(B), 1))
    d = np.zeros((len(B), 1))
    X = np.zeros((len(B), 1))
    c[0] = A[0,1] / A[0,0]
    d[0] = B[0] / A[0,0]
    for i in range(1, len(B)-1):
        c[i] = A[i,2] / (- A[i,0] * c[i-1] + A[i,1])
        d[i] = (B[i] - A[i,0]*d[i-1]) / (- A[i,0] * c[i-1] + A[i,1])
    X[-1] = (B[-1] - A[-1, -2]*d[-2]) / (A[-1,-1] -A[-1,-2]*c[-2] )
    for i in range(len(B)-2, -1, -1):
        X[i] = - c[i]*X[i+1] + d[i]
    return X
