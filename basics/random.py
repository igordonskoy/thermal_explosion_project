def Fk_stuffe(Fk_avg, N, delta, position):
    if position < 0.5:
        p1 = position
        p2 = position + delta
        if p2 > 1:
            print('Positioning is poorly chosen.')
    elif position >= 0.5:
        p1 = position - delta
        p2 = position
        if p1 < 0:
            print('Positioning is poorly chosen.')
    else:
        print('Positioning is poorly chosen.')
    i1 = round(p1*N)
    i2 = round(p2*N)
    Fk_max = Fk_avg/delta
    Fk = np.zeros((N,1))
    Fk[i1:i2] = Fk_max
    return Fk


def Fk_distr(Fk_avg, N, sigma=0, mode=0):
    if mode==0:
      noise = np.random.randn(N, 1)               # this code is for
      Fk = Fk_avg + sigma*noise                   # N-distribution
      Fk[Fk < 0] = 0                              # but with
      Fk[Fk > 2*Fk_avg] = 2*Fk_avg                # cutted tails
    if mode == 1:
      noise = 2*(np.random.rand(N,1) - 1/2)
      Fk = Fk_avg + sigma * noise
    if mode == 2:
      noise = np.random.randint(2, size=(N,1))   # the code for
      Fk_i = Fk_avg * N / noise.sum()            # 0/1 distribution
      Fk = noise * Fk_i                          # with normalization
    return Fk


def random_Y(Y_avg, NX, NY):  # bins filling
    dY = Y_avg/NY
    Y_random = np.zeros((NX,1))
    k = 0
    while k < NY*NX:
        k += 1
        i = np.random.choice(range(NX))
        Y_random[i] += dY
    S = np.trapz(Y_random, dx=1/(NX-1), axis=0)
    dS = Y_avg - S
    Y_random += dS
    return Y_random


