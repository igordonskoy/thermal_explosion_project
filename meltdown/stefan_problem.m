clear

% == Stefan problem in standard bodies ==
% == heated by outer boundary ==

n = 1;                                                                      % symmetry [0 - flat, 1 - cylinder, 2 - sphere]
L = 1; N = 100; h = L/N; x = h:h:L;                                         % spatial grid
[D,Cl,Cr] = transfer_matrices(N);                                           % 
Dn = diag([0 x(1:end-1)].^n)*Cl + diag(x.^n)*Cr; En = diag(x.^n);           % curvature
T = repmat(-1, N, 1);                                                       % initial condition (temperature distribution)
kappa = 1; c = 1.5;                                                         % ratio of phases conductivities
Tenv = 2; Tm = 0; St = 10; Bi = 10;                                         % heating and melting temperature; Stefan dn Biot numbers
s = L; im = N;                                                              % melting front initial condition
Y = [T; ];                                                                  % recording
dt = 0.1e-2; t = 0; t_end = 10; melting_period = 0; k = 1;                  % time grip
tic
while t(end) < t_end                                                        % time cycle
    k = k + 1; t(k) = t(k-1) + dt;                                          % time update
    qm = 0;                                                                 % heat flux initialization
    B = En*T;                                                               % previous temperature
    if melting_period                                                       % if sample melts
        if im > 1 & im < N                                                  % and we are not an boundaries
            aT = diag([repmat(1,(im-1),1); repmat(kappa,(N-im+1),1)]);      % then we construct
            A = En - dt/h^2*(aT*Dn);                                        % transition matrix
            A(N,N) = 1 + h*Bi; A(N,N-1) = -1;  B(N) = h*Bi*Tenv;            % with outer boundary condition 
            A(1,1) = -1; A(1,2) = 1;B(1) = 0;                               % and inner boundary condition
            A(im,im) = 1; A(im,im-1) = 0; A(im,im+1) = 0; B(im) = Tm;       % and Stefan boundary condition
            T = progonka(A,B);                                              % then we solve it
            qm = (T(im+1)-Tm)/h - (Tm - T(im-1))/h;                         % and calculate heat fluxes at phase boundary
        elseif im == N                                                      % if we are at the outer (heated) boundary
            A = En - dt/h^2*Dn;                                             % then transfer matrix is similar
            A(N,N) = 1; B(N) = Tm; B(1) = 0;                                % but the boundary temperature is the melting temperature
            T = progonka(A,B);                                              % and solve it
            qm = Bi*(Tenv-Tm) - (Tm - T(im-1))/h;                           % again, calculate heat flux difference
        elseif im == 1                                                      % finally, if the boundary is close to centre
            A = En - dt/h^2*(kappa*Dn);                                     % we construct transition matrix full of liquid conductivity
            A(N,N) = 1 + h*Bi; A(N,N-1) = -1; A(1,1) = 1;                   % then outer boundary condition is as usual
            B(N) = h*Bi*Tenv; B(1) = Tm;                                    % but inner temperature is the melting temperature
            T = progonka(A,B);                                              % the solve it
            qm = (T(im+1) - Tm)/h;                                          % anyway, calculate heat flux difference
        end
    else                                                                    % what if sample doens not melt yet/already?
        A = En - dt/h^2*(1*(im == N) + kappa*(im == 1))*(Dn);               % we solve heat transfer problem
        A(N,N) = (1*(im == N) + kappa*(im == 1)) + h*Bi;                    % just checking
        A(N,N-1) = -(1*(im == N) + kappa*(im == 1));                        % if it yet or already
        A(1,1) = -1; A(1,2) = 1; B(N) = h*Bi*Tenv; B(1) = 0;                % without phase boundary condition
        T = progonka(A,B);                                                  % solve it
    end
    uf(k) = qm/St;                                                          % heat flux difference determines phase front velocity
    s(k) = max([s(k-1) - dt*uf(k), 0]);                                     % updating front location
    if ~melting_period & s(k) & max(T) >= Tm                                % if melting has not started,
        melting_period = 1;                                                 % it may start
    end
    if melting_period & ~s(k)                                               % if sample completely melted 
        melting_period = 0;                                                 % then 
        t_end = t(end);                                                     % stop calculation
    end
    
    [xm,im] = near(s(k), x);                                                % allocate the melting front
    Y(:,k) = T;                                                             % record temperature
end
toc   
% == Plotting ==
figure
subplot(1,2,1)
plot(x,T,x(im),Tm,'o')
subplot(1,2,2)
plot(t,s)
