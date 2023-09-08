import numpy as np
import numpy.matlib
from matplotlib import pyplot as plt

def film_ignition_Q(N, Fk, Le, s, Ys, Ym, Ar):
    h = 1/(N-1); x = np.linspace(0,1,N);
    yeq = np.exp(s*Ys); # saturation (Ys should be less than 0)
    Pe = 0;
    if B > 0:
        Pe = 1/Le*np.log(1 + yeq/(1-yeq)); # Stefan flow
    else:
        Pe = 1/Le*yeq; # straight diffusion flow
    y = yeq/(np.exp(Pe) - 1)*(np.exp(Pe) - np.exp(Pe*x));
    Y0 = np.matlib.repmat(Ys,N,1);
    A = np.zeros((N,N)); #heat flow matrix
    A[0,0] = 1;
    for i in range(1,N-1):
        A[i,i-1] = -1 - h*Pe; A[i,i] = 2 + h*Pe; A[i,i+1] = -1;
    A[N-1,N-1] = 1;
    k = 0;
    dF = 1; Q = False;
    while dF > 1e-3:
        k = k + 1; S = np.zeros((N,1));
        for i in range(0,N):
            S[i] = h**2*Fk*y[i]*np.exp(Y0[i]/(1 + Ar*Y0[i]));
        #S = h**2*Fk*np.multiply(y,np.exp(Y0/(1 + Ar*Y0))); # reaction heat
        #print(S)
        B = S; B[0] = Ys; B[N-1] = Ym;
        Y = np.linalg.solve(A,B); #solving
        dF = max(abs(Y - Y0)); #checking convergence
        if k > 100:
            dF = 0; #too long - apparently, unstable
        if max(Y) > 1e3:
            Q = True; #explosion
        Y0 = Y;
    return Q # bool
