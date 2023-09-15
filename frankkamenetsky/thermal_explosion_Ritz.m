function [B,Y1,a] = thermal_explosion_Ritz(Fk,Ar,Bi,Np)
c1 = 1; c2 = 1; c3 = 1;
for i = 1:(Np+1)
    for j = 1:(Np+1)
        D0(i,j) = 0;
        if (i+j-3) ~= 0
            D0(i,j) = (i-1)*(j-1)/(i+j-3);
        end
    end
end
b = [0:Np];
a = zeros(Np+1,1);

dx = 0.001;
x = stolb(0:dx:1);
Y = zeros(length(x),1);
for i = 1:Np
    Y = Y + a(i)*x.^(i-1);
end
B = 0;
dY = 10; p = 0;
while dY > 1e-3
    p = p + 1;
    Y0 = Y;
    q = Fk*exp(Y0./(1+Ar*Y0));
    for k = 1:(Np+1)
        s = q.*x.^(k-1);
        S(k,1) = trapz(x,s);
    end
    D = [ D0 (b+Bi)' ; (b+Bi) 0];
    D(end+1,2) = 1; D(1:Np+1,end+1) = 1;
    w = D\[S; 0;0];
    a = w(1:Np+1);
    Y = zeros(length(x),1);
    for i = 1:(Np+1)
        Y = Y + a(i)*x.^(i-1);
    end
    dY = max(abs(Y-Y0));
    if p > 1000 | max(Y) > 1e1
        B = 1; Y1 = Y;
        break
    end
end
toc
Y1 = Y;
