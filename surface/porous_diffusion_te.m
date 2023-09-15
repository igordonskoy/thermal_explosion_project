clear
% The program for solving the problem of heat and mass transfer in porous
% media with outer flow of combustible gas 

N = 500;                                % grid points number
L = 1; h = L/(N-1);                     % grid step
x = 0:h:L;                              % spatial coordinate

[A] = transfer_matrices(N);             % 1D transfer matrix
bc = 1;                                 % boundary condition for x = 1

E = eye(N);                             % identity matrix

T = zeros(N,1); C = ones(N,1);          % initial conditions
T0 = 0; C0 = 1;                         % boundary values

t = 0; dt = 2e-3; t_end = 10;           % time grid
T_rec = T; C_rec = C;                   % recording variables
Da0 = 2e-2; Q = 50; Ar = 2e-2; Le = 0.5;  % Damkohler number, adiabatic temperature, Arrhenius number, Lewis number
Da_margin = 0;                          % Damkohler number relative margin (from 0 to 1)
Da = Da0*(1 + Da_margin*(2*rand(N,1) - 1));   % random distribution with average <Da> = Da0
if bc == 3                              
    Bi = 100;                           % Biot number (if boundary conditions are of the third kind)
end
ST = E - dt/h^2*A;                      % heat transfer transition matrix
SC0 = E - 1/Le*dt/h^2*A;                % diffusion transition matrix (without chemical reaction)
ST(1,1) = 1; ST(1,2) = -1;              % boundary condition of the second kind for me 
SC0(1,1) = 1; SC0(1,2) = -1;            % and for my son too
if bc == 3                              % boundary condition of the third kind
    SC0(N,N) = 1 + h*Bi; SC0(N,N-1) = -1; 
    ST(N,N) = 1 + h*Bi; ST(N,N-1) = -1;
end
if bc == 1                              % boundary condition of the first kind
    SC0(N,N) = 1; ST(N,N) = 1;
end
k = 1;                                  % counter
maxTrate = -1; t_ign = -1; T_ign = -1;  % variable to store maximum heating rate, ignition time and ignition temperature
while t(end) < t_end
    k = k + 1; t(k) = t(k-1) + dt;      % updating time
    C_prev = C; T_prev = T;             % previous solution
    K = diag(Da.*exp(T./(1 + Ar*T)));    % calculating reactivities
    K(1) = 0; K(end) = 0;               % excluding boundaries
    SC = SC0 + dt*K;                    % diffusion transition matrix (with chemical reaction)
    BC = C_prev; BC(1) = 0; 
    BT = T_prev; BT(1) = 0; 
    if bc == 3                          % boundary condition of the third kind
        BC(end) = h*Bi*C0; BT(end) = h*Bi*T0;
    end
    if bc == 1                          % boundary condition of the first kind
        BC(end) = C0; BT(end) = T0;
    end
    C = progonka(SC, BC);               % solving diffusion problem under fixed T
    dC = dt*K*C;                        % calculating chemical sources
    dT = dC*Q; dT(1) = 0; dT(end) = 0;  % calculating temperature jumps
    BT = BT + dT;                       % adding them to temperature
    T = progonka(ST, BT);               % solving heat transfer problem under fixed C
    
    T_rec(:,k) = T; C_rec(:,k) = C;     % recording
    if k > 1
        if (T_rec(1,k) - T_rec(1,k-1))/dt > maxTrate        % cheking for ignition
            maxTrate = (T_rec(1,k) - T_rec(1,k-1))/dt;      % specifically, maximum heating rate point
            t_ign = t(k); T_ign = T_rec(1,k);               % recording current extremum
        end
    end
    if any(T<0) | any(C<0)              % checking for negativity
        disp('Unphysical solution.')
        break
    end
    ax = 0;
    
end
    
figure
if t_ign == -1 | T_ign < 1
    plot(t, T_rec(1,:))
    title('Center temperature dynamics: no ignition')
else
    plot(t, T_rec(1,:), t_ign, T_ign, 'o', [t_ign, t_ign], [T_ign, T0], '--')
    title(['Center temperature dynamics: \tau_{ign} = ', num2str(t_ign)])
end
xlabel('Time')
ylabel('\theta(0)')

figure
subplot(1,2,1)
plot(x, T)
title('Final temperature distribution')
xlabel('Spatial coordinate')
ylabel('\theta')
subplot(1,2,2)
plot(x, C)
title('Final concentration distribution')
xlabel('Spatial coordinate')
ylabel('\phi')

