clear

N = 200; L = 0.5;
h = L/(N-1);
x = 0:h:L;

A = zeros(2*N,2*N);
B = zeros(2*N);

Fk1 = 1.28; Fk2 = Fk1;
Pe1 = 1.0; Pe2 = Pe1;
a = 1;
Bi = 0.1;

Y0 = zeros(N,1); Z0 = zeros(N,1);
T0 = [Y0; Z0];

A(1,1) = 1; A(N,N-1) = 1; A(N,N) = -2; A(N,2*N) = 1;
A(N+1,N+1) = -1; A(N+1,N+2) = 1; 
A(2*N,2*N) = -2; A(2*N,N) = 1; A(2*N,2*N-1) = 1;
% A(2*N,2*N) = 1; A(2*N,N) = -1;
for i = 2:(N-1)
    A(i,i-1) = 1 + h*Pe1; A(i,i+1) = 1;
    A(i,i) = -2 - h^2*Bi - h*Pe1; A(i,N+i) = h^2*Bi;
    A(N+i,N+i-1) = a; A(N+i,N+i+1) = a + h*Pe2;
    A(N+i,N+i) = -2*a - h^2*Bi - h*Pe2; A(N+i,i) = h^2*Bi;
end
is_stat = 0;
while is_stat == 0
    Sy = Fk1*exp(Y0); Sz = Fk2*exp(Z0);
    B = -h^2*[Sy; Sz];
    B(1) = 0; B(N) = 0; B(N+1) = 0; B(2*N) = 0;
    T = A\B;
    max_diff = max(abs(T-T0));
    T0 = T;
    Y0 = T(1:N); Z0 = T(N+1:2*N);
    q = -Bi*(Y0 - Z0);
    if max_diff < 1e-3
        is_stat = 1;
        Q = 0; fig_title = 'No babakh';
    end
    if max(T) > 1e3 | sum(isnan(T)) ~= 0 
        Q = 1; fig_title = 'Babakh';
        break
    end
    disp(num2str(max(T)))
end
figure
plot(x,Y0,x,Z0)
title(fig_title)
if Bi ~= 0
    figure
    plot(x,q)
    title('heat flow')
end
