def transfer_matrix_col(N, lbc=1):
    A = np.zeros((N, 3))
    E = np.zeros((N, 3))
    E[0,0], E[-1,-1] = 1, 1
    for i in range(1,N-1):
        A[i,1] = -2
        A[i,0] = 1
        A[i,2] = 1
        E[i,1] = 1
    A[0,0] = 1
    if lbc == 2:
        A[0,1] = -1
    A[N-1,-1] = 1
    if lbc == 0:
        A[0,0] = 0
        A[N-1, -1] = 0
    
    return A, E
