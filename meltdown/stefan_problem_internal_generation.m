clear
% == Program for modelling Stefan problem =====
% == in a body with internal heat generation ==
% == Based on [Crepeau et al., 2008] ==
% =========================================================================
% == Parameters ==
n = 1;                                                                      % geometry [0 - flat, 1 - cylinder, 2 - sphere]
N = 200;                                                                    % the number of spatial nodes
Tin = 0; Tm = 1; Tenv = 0;                                                  % characteristic temperatures [initial, melting, environment]
kappa = 2; W = 4;                                                           % phase conductivities ratio and heat generation power   
Bi = 10; St = 10; Bo = 0; T_0 = 1;                                           % Biot, Stefan and Boltzmann numbers; ground zero radiation level, T0/(Tm - T0)
% == Problem initialization ==
L = 1; h = L/N; x = h:h:L;                                                  % spatial grid
T = repmat(Tin, N, 1);                                                      % initial distribution of temperature
X = ones(N, 1);                                                             % initial distribution of solid phase
E = [[1, 0, 0]; repmat([0, 1, 0], N-2,1); [0, 0, 1]];                       % identity matrix
t_end = 20; t = 0; k = 1; s = 0;                                            % full time and starting condition
dt = 1e-3; Fo = dt/h^2;                                                     % time step and cell Fourier number
% == Time cycle ==
RES1(:,1) = [T; X];                                                         % track solution
while t(end) < t_end                                                        % stop condition
    k = k + 1; t(k) = t(k-1) + dt;                                          % new time node         
    lambda = 1*X + kappa*(1-X);                                             % conductivities
    lambda = (lambda + [0; lambda(2:end)])/2;                               % average conductivities
    A = zeros(N,3);                                                         % transfer matrix
    for i = 2:(N-1)
        phi = (x(i-1)/x(i))^n;                                              % geometric factor
        A(i,:) = [-phi*lambda(i), (1*lambda(i+1) + phi*lambda(i)), -1*lambda(i+1)];
    end
    T0 = T; X0 = X;                                                         % double previous moment solution
    Q = E + Fo*A; B = T0 + dt*W;                                            % transition matrix and initial condition
    Q(1,1) = 1 + Fo*lambda(2); Q(1,2) = -Fo*lambda(2);                      % inner boundary condition
    Tr = (1 + T0(N)/T_0);                                                   % radiation temperature
    Q(N,2) = -1*lambda(N);                                                  % outer boundary condition
    Q(N,3) = 1*lambda(N) + h*Bi + 4*h*Bo*Tr^3;                              % conductive, convective and radiative heat flux
    B(N) = h*Bi*Tenv + h*Bo*(1 - Tr^4 + 4*Tr^3*T0(N));                      % 
    T = progonka_col(Q,B);                                                  % solve heat transfer problem as if there were no phase transitions
    q_out(k) = Bi*(T(N) - Tenv) + Bo*((1 + T(N)/T_0)^4 - 1);
    for i = 1:N
        if (T(i) > Tm) & (X0(i) > 0)                                        % melting condition
            deltaT = T(i) - Tm;                                             % overheating
            deltaX = deltaT/St; X(i) = X0(i) - deltaX; T(i) = Tm;           % melting by heat balance
            if deltaX > X0(i)                                               % preventing overmelting 
                T(i) = Tm + (deltaT - St*X0(i));  X(i) = 0;                 % melting by mass balance
            end
        end
        if (T(i) < Tm) & (X0(i) < 1)                                        % crystallization condition
            deltaT = T(i) - Tm;                                             % overcooling
            deltaX = -deltaT/St; X(i) = X0(i) + deltaX; T(i) = Tm;          % crystallization by heat balance
            if deltaX > (1 - X0(i))                                         % preventing overcrystalization
                T(i) = Tm - deltaT + St*X0(i); X(i) = 1;                    % crystallization by mass balance
            end
        end
    end
    s(k) = s(k-1);                                                          % phase transition front coordinate (discrete)
    if any(X == 0)                                                          % if there is melt
        xl = x(X > 0);                                                      % then its terminate coordinate
        s(k) = xl(1);                                                       % is the phase surface
    end
    
    if s(k) >= 1                                                            % oh no, we have a meltdown!
        disp('Meltdown')                                                    % signalizing
        t_end = t(k);                                                       % stop calculation
    end
%     RES1(:,k) = [T; X];                                           %     recording
end

% == Plotting ==

figure
plot(x,T, x, X)
xlabel('Spatial coordinate')
legend({'Temperature', 'Solid fraction'})

figure
plot(t, s)
xlabel('Time')
ylabel('Phase transition front location')

% figure
% subplot(1,2,1)
% plot(x, RES1(1:N,:), [x(1), x(end)], [Tm, Tm],'--')
% subplot(1,2,2)
% plot(x, RES1(N+1:2*N,:))
