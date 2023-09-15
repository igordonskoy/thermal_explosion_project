function [q, Y] = thermal_explosion_dE(Fk,Ar,dE,N,Bi,Pe);

dz = 1/(N-1);
z = 0:dz:1;
Y = zeros(N,1);

st = 4*dE; ds = st/50;
s = [-st:ds:st];
g = exp(-s.^2/2/dE^2);
C = trapz(s,g); f = g/C;
A = zeros(N,N); B = zeros(N,1);
A(1,1) = 1; A(1,2) = -1; B(1) = 0;
A(N,N-1) = 0; A(N,N) = 1; B(N) = 0;
for i = 2:(N-1);
    A(i,i-1) = 1 ; A(i,i+1) = 1 + dz*Pe;
    A(i,i) = -2 - dz*Pe - dz^2*Bi;
end

v = 1; k = 0;
while v > 0
    k = k + 1;
    Y0 = Y;
    c3 = zeros(N,1);
    F = 1;
    for i = 2:(N-1)
        fq = f.*exp(-s/Ar).*exp(s*Y0(i)/(1+Ar*Y0(i)));
        F = trapz(s,fq);
        Z = Y0(i)/(1+Ar*Y0(i));
        S(i)= Fk*exp(Z)*F;
    end
 
    B(2:end-1) = -dz^2*S(2:end);
    Ys = progonka(A,B);
    ax = 0;
    if k > 1
        dYmax = max(abs(Ys - Y0));
        if dYmax < 1e-3
            v = 0;
            q = 1;
            disp(['No babakh'])
        end
    end
    Y = Ys;
    if k > 1e6
        disp('Babakh!')
        disp(['dYmax = ',num2str(dYmax)])
        v = 0;
        q = 0;
    end
    if sum(isnan(Ys)) ~= 0
        disp('Babakh!')
        v = 0;
        q = 0;
    end
    if max(Ys) > 1e1
        disp('Babakh!')
        v = 0;
        q = 0;
    end
end
