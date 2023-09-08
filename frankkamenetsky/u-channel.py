def regenerative_ignition_Q(Fk, N=100, L=1, Pe=0, Bi=0, vg=False):
    h = L/(N-1)
    x = np.linspace(0,L,N)
    A = np.zeros((N, N))
    B = np.zeros((N,1))
    Y0 = np.zeros((N,1))
    A[0,0] = -1 - h*Pe
    A[0,1] = 1
    A[-1,-1] = -1
    A[-1,-2] = 1
    for i in range(1,N-1):
        A[i,i-1] = 1 + h*Pe
        A[i,i+1] = 1
        A[i,i] = -2 - h*Pe - 1*h*h*Bi
        A[i,-i-1] += h*h*Bi
    is_solved = False
    while not is_solved:
        S = Fk*np.exp(Y0)
        B[1:N-1] = -h*h*S[1:N-1]
        Y = np.linalg.solve(A, B)
        dY = max(abs(Y-Y0))
        if dY < 1e-3:
            is_solved = True
            Q = False    # no ignition
        if max(Y) > 1e2:
            is_solved = True
            Q = True    # ignition detected
        Y0 = Y
    if vg:
        fig, a1 = plt.subplots()
        a1.plot(x, Y)
        a1.grid()
    return [Q, x, Y0]
