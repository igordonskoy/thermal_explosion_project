clear

Le = 1; Tm = 10; s = 10;
X = 0:0.01:10;
Ar = 1e-2;
Y1 = 1e-3; Y2 = 0.05; N = 1000; dy = (lg(Y2)-lg(Y1))/N; y = lg(Y1):dy:lg(Y2);
Y = 10.^[y];
Y = 0.042;
Q = zeros(length(Y),3);

for j = 1:length(Y);
    Da = Y(j);
    Keq = exp(s*X);
    c = (1 + Da*exp(X*(1-s)))./(1 + Da*exp(X).*(1 + exp(-s*X)));
    A = (s/Ar./(1+Ar*X) + ln(c./(1-c)));
    w = Da*exp(X).*(1 + exp(-s*X)).*((1+Da*exp(X*(1-s)))./(1 + Da*exp(X).*(1 + exp(-s*X))) - 1./(1+exp(s*X)));
    S = X.^2./(1+Ar*X) + 1/Le*(1-c).*ln(c) + 1/Le*(s./(1+Ar*X) + Ar*ln(c./(1-c))).*w;
    f = Tm/Le*Da*exp(X).*(1 + exp(-s*X)).*((1+Da*exp(X*(1-s)))./(1 + Da*exp(X).*(1 + exp(-s*X))) - 1./(1+exp(s*X)));
    F = X - f;
    k = 0;
    for i = 2:length(X)
        if F(i-1)*F(i) < 0
            k = k + 1;
            Q(j,k) = (X(i)+X(i-1))/2;
            u(k) = i;
        end
    end
end
Q = [Y' Q];
figure
plot(X,X-f,[X(1) X(end)],[0 0],X(u),zeros(1,length(u)),'o')

figure
plot(X,S,X(u),S(u),'o')
