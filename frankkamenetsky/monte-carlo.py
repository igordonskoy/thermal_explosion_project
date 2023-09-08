def find_critical_Fk_stuffe(delta, position, N, max_iter):
    Fk1 = 0
    expl_Fk2 = False
    Fk2 = 0.5
    while expl_Fk2==False:
        Fk2 *= 2
        Fk_p = Fk_stuffe(Fk2, N=N, delta=delta, position=position)
        _, expl_Fk2 = given_Fk_solve(Fk_p, N=N, max_iter=max_iter, view_plot=False)
    k = 0
    Fks = (Fk1 + Fk2)/2
    while (Fk2-Fk1) > 1e-4*Fk2:
        k += 1
        Fks = (Fk1 + Fk2)/2
        Fk_p = Fk_stuffe(Fks, N=N, delta=delta, position=position)
        _, expl_Fks = given_Fk_solve(Fk_p, N=N, max_iter=max_iter, view_plot=False)
        if expl_Fks:
            Fk2 = Fks
        else:
            Fk1 = Fks
        if k > max_iter:
            print('Iteration limit is exceeded.')
            break
    return Fks


def given_Fk_solve(Fk_Y, N=100, max_iter=1000, view_plot=False):
    T = np.zeros((N, 1))    # temperature distribution
    B = np.zeros((N,1))
    dQ = 1e6
    expl = False
    epsQ = 1e-4
    k = 0
    h = 1/(N-1)
    A = transfer_matrix(N, lbc=2)
    while dQ > epsQ:
        T0 = T
        S = h**2 * Fk_Y* np.exp(T0)
        B[2:N-1] = - S[2:N-1]
        T = progonka(A, B)
        dQ = max(abs(T - T0))
        k += 1
        if k > max_iter:
            dQ = 0
        if max(T) > 10:
            dQ = 0
            expl = True
    if view_plot:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        fig.tight_layout(pad=5.0)
        ax1.plot(np.linspace(0, 1, N), Fk_Y)
        ax1.set_xlabel('x')
        ax1.set_ylabel('Fk')
        ax1.set_title(f'mean(Fk) = {round(Fk_Y.mean(),2)}')
        ax2.plot(np.linspace(0, 1, N), T)
        ax2.set_title(f'max(T) = {round(T.max(),3)}, iter = {k}')
        ax2.set_xlabel('x')
        ax2.set_ylabel('T')
        plt.show()
    return T, expl


def random_Fk_MC(Fk_avg, N_MC=10, Nx=500, max_iter=100, sigma=0):
    explosions = np.zeros((N_MC, 1))
    T_max = np.zeros((N_MC, 1))
    for i in range(N_MC):
        #Fk_given = random_Y(Fk_avg, Nx, 100)
        Fk_given = Fk_distr(Fk_avg, Nx, sigma=sigma)
        T_st, explosion_indicator = given_Fk_solve(Fk_given, N=Nx, max_iter=max_iter, view_plot=False)
        T_max[i] = T_st.max()
        explosions[i] = explosion_indicator
    return explosions, T_max



