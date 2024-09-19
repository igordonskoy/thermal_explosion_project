% Program to calculate nonstatinary ignition of quescent cylindrical gas volume
% around a coil heated with the constant current

clear

% = Spatial grid =
r0 = 1; D = 2;              % Coil radius is 1, D is a radius of gas volume
R = 1 + D;                  % Overall radius
h = 1e-3; h2 = h^2;       % Sparial grid step
r = 0:h:R; N = length(r);   % grid and its size
Ns = round(N*r0/R);         % surface element
phi = r ./ (r + h);         % cylindrical ratio

% = Thermophysical constants (numbers) =
Bo = 1; Ar = 1e-2; Le = 1;  % Boltzmann number (for coil surface), Arrhenius number, Lewis number
Da0 = 1e-2; Tad = 50;       % Damkohler number (for ground temperature), adiabatic temperature
T_out = 0; Y_out = 1;       % environment temperature (default zero), environment reagent fraction (default 1)
gamma = 100; kappa = 100;   % heat capacity and conductivity of coil material (ratio to gas properties)
Q = 2;                      % Joule heat generation intensity

% = Time grid =
t_end = 5;                  % full time (in relative units, i.e. Fourier number)
dt0 = 0.5e-3;                 % time step
N_rec = 200; dt_rec = t_end/N_rec; t_rec = 0; % recording parameters
ignition = false;           % ignition check

% = There we go for matrices =
Fo_w = kappa/gamma*1/h2;    % coil cell Fourier number divided by dt
Fo_g = 1/h2;                % gas cell Fourier number divided by dt
q = zeros(N,1);             % heat sources distribution
[A, ET] = transfer_matrices_col(N); [A, ED] = transfer_matrices_col(N-Ns);  % matrices are ready
A = zeros(N,3); B = zeros(N,1); % this line is just to fix dimensions
ET(1,1) = 1; ET(1,2) = -1;  % Neumann conditions at the coil center
ET(N,3) = 1; %              % Dirichlet conditions at the outer boundary
for i = 2:(Ns-1)            % filling the transfer matrix inside the coil
    A(i,1) = -Fo_w*phi(i); A(i,3) = -Fo_w;
    A(i,2) = (1+phi(i))*Fo_w;
end
A(Ns,1) = kappa; A(Ns,3) = 1;   % surface conditions: difference of heat fluxes
A(Ns,2) = -(1+kappa); B(Ns) = 0; ET(Ns,2) = 0;  % zero for identity matrix
C = zeros((N-Ns),3);        % diffusion matrix out of coil only
ED(1,1) = 1; ED(1,2) = -1; ED(end,3) = 1;   % conditions are the same as for temperature
for i = (Ns+1):(N-1)        % filling the transfer matrices outside the coil
    A(i,:) = [-phi(i)*Fo_g, (1+phi(i))*Fo_g, -Fo_g];                % heat transfer
    C(i-Ns+1,:) = [-phi(i)*Fo_g/Le, (1+phi(i))*Fo_g/Le, -Fo_g/Le];  % diffusion
end
C(N-Ns,:) = [0, 0, 0];      % the last line should be empty (didn't figure out why isn't it already zero)

% = Initial conditions, surface values =
T = zeros(N,1); Ts = T(Ns);
Y = ones(N-Ns,1); Ys = Y(1);   

% = Preparing to start = 
t = 0; t_rec = 0; dt = dt0; % zeros and input time step
k = 1; VAR.T(:,k) = T; VAR.Y(:,k) = Y;  % counting and recording initial conditions
A0 = A; C0 = C;             % saving untouched matrices (we're gonna touch them latter)
% = Main cycle =
while t(end) < t_end        % checking time
    k = k + 1;              % refresh count variable
    T0 = T; Y0 = Y;         % saving previous results
    Da = Da0*exp(T(Ns+1:N)./(1+Ar*T(Ns+1:N)));  % kinetic coefficients in the gas area
    Dam = min([max(Da)*Y(Ns), 1e2]);                  % finding out the biggest one (with the cutting though)
    if k > 3                        % a simple way to adapt time step
%         dt = dt0/(1 + Dam*Y(Ns));   % it should switch off when the gas already ignited
        dt = dt0/(1 + Dam);   % it should switch off when the gas already ignited
    end
    % == Splitting consequence ==
    % === Diffusion subproblem (T(r) is fixed) ===
    C = C0; C(2:end-1,2) = C(2:end-1,2) + Da(2:end-1);  % linear kinetics is nice, we can easily include it into diffuston matrix
    G = Y0; G(1) = 0; G(end) = Y_out;   % zero gradient right, main flow left
    QD = ED + dt*C;                     % transfer matrix as it is
    Y = progonka_col(QD, G);            % TDMA solution
    w = Da.*Y;                          % evaluating chemical reaction rates 
    % === Heat subproblem (C(r) is fixed) ===
    q(1:Ns-1) = Q;                      % Joule heating in coil
    q(Ns+1:N) = Tad*w;                  % exothermic heating in gas
    q_rad = Bo*((1+Ar*T(Ns))^4 - 1);    % radiative heat flux at the surface 
    B = T + dt*q; B(Ns) = h*q_rad; B(1) = 0; B(N) = T_out;  % zero gradient right, main flow left
    A = A0;     % we're gonna change one line
    A(Ns,2) = A(Ns,2) - h*4*Bo*Ar*(1+Ar*Ts(k-1))^3; % surface radiative flux linearized
    B(Ns) = B(Ns) - h*4*Bo*Ar*(1+Ar*Ts(k-1))^3*Ts(k-1);             % its second part
    QT = ET + dt*A;                     % heat transfer matrix as it is 
    T = progonka_col(QT,B);             % TDMA again
    q_out = (T(N-1) - T(N))/h;          % conductive heat flux into outer space
    % == Processing the solution ==
    Ts(k) = T(Ns); Ys(k) = Y(1);        % surface values
    qw1 = kappa*(T(Ns-1)-T(Ns))/h; qw2 = (T(Ns)-T(Ns+1))/h; % surface heat fluxes
    q_b(k,:) = [q_rad, q_out, qw1, qw2];    % recording surface heat balance
    t(k) = t(k-1) + dt;                 % updating time
    if ~ignition                        % if there wasn't ignition
        if min(Y) < 1e-3                % then check it: the criterion is reagent depletion
            ax = 0;                     % line to stop in debugging mode
            ignition = true;            % if reagent is depleted, then ignition is succesfull
            t_ign = t(k);               % record ignition time
            Q = 0;                      % turn off the heating
            Y_out = 0;                  % turn off gas supply
            VAR.T(:,end+1) = T; VAR.Y(:,end+1) = Y; % record temperature and concentration profiles
            t_rec(end+1) = t(k);        % record time stamp
        end
    end
    if t(k) >= t_rec(end) + dt_rec;     % recording profiles in any case
        VAR.T(:,end+1) = T; VAR.Y(:,end+1) = Y; % profiles
        t_rec(end+1) = t(k);            % time stamps
    end
end
% = Plotting =
figure  % last temperature profile
plot(r, T, r(Ns+1:N), Y, [r(Ns), r(Ns)], [0, max([max(T), max(Y)])], '--k')
title('Final distribution')
xlabel('r')
ylabel('T, Y')

if ignition
    figure % ignition shot
    plot(r,VAR.T(:,t_rec == t_ign), r(Ns+1:N), VAR.Y(:,t_rec == t_ign), [r(Ns), r(Ns)], [0, Tad], '--k') 
    title(['Ignition moment: ', num2str(t_ign)])
    xlabel('r')
    ylabel('T, Y')
end

figure % phase trajectory piece for surface
plot(Ts, Ys, '-', Ts(t == t_ign), Ys(t == t_ign), 'o')                           
title('Surface dynamics')
xlabel('Temperature')
ylabel('Concentration')
