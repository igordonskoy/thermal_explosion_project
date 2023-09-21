% Program for solving unstatioanry Frank-Kamenetskii problem with
% stochastic boundary temperature
clear

N = 300;                                    % the number of grid points
h = 1/(N-1);                                % grid step

Fk = 0.7; Ar = 1e-2;                        % thermochemical parameters
T_env = 0; Bi = 100;                        % boundary parameters
sigma = 0.1;                                % variance
t_end = 100; dt = 0.001;                    % time grid

E = eye(N); E(1) = 0; E(N) = 0;             % choped identity 
T_in = zeros(N,1); T = T_in;                % initial temperature distribution
T_c = T(1); T_b = T(N); T_m = max(T);       % recording
A = transfer_matrices(N);                   % tdm
A(1,1) = 0; A(1,2) = 0;                     % nullification of boundaries
A(N,N) = 0; A(N,N-1) = 0;                   %
t = 0; k = 1;                               % running variables
% T_rec = T;                                  % recording dump
solved = false;                             % solution flag
while ~solved 
    k = k + 1; t(k) = t(k-1) + dt; T0 = T;  % counting and refreshing
    T_env(k) = T_env(k-1) + sqrt(dt)*sigma*randn;      % temperature drift
    Q = E - dt/h^2*A;                       % transition matrix
    Q(1,1) = 1; Q(1,2) = -1;                % zero gradient at center
    Q(N,N) = 1 + h*Bi; Q(N,N-1) = -1;       % third kind bc
    S = dt*Fk*exp(T0./(1 + Ar*T0));         % chemical sources
    B = T0 + S; B(1) = 0; B(N) = h*Bi*T_env(k); % temperature update
    T = progonka(Q,B); % T_rec(:,k) = T;    % solution/recording
    if max(T) > 10 | t(end) >= t_end        % either explosion or time limit
        solved = true;                      % ends the calculations
    end
    T_c(k) = T(1); T_b(k) = T(N); T_m(k) = max(T);  % record characteristic temperatures
end
% Plotting
figure
plot(0:h:1, T)                                  % final temperature distribution along the sample
xlabel('\xi')
ylabel('\theta')

figure                                          % dynamics of characteristic temperatures
plot(t, T_c, t, T_b, t, T_env, 'x', t, T_m, '.')
title(['\tau_{end} = ', num2str(t(end))])
xlabel('\tau')
ylabel('\theta')
legend({'\theta_{center}', '\theta_{bound}', '\theta_{env}', '\theta_{max}'})
