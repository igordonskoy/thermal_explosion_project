clear
% Thermal explostion in Semenov reactor with DAEM kinetics

isothermal = false;                         % flag to exclude/include temperature equiation

sD = 0.04; s = -0.99:0.001:0.99;            % variance and phase space grid
f1 = exp(-s.^2/2/sD^2); Cf = sum(f1);       % PDF for DAEM
f = 1/Cf*f1; g = f;                         % normalization
Ar = 1e-2; Tm = 10;                         % Arrhenius factor and adiabatic temperature

N = f; NS = sum(f); Nb = 1;                 % initialization of pseudocomponents
T = 0; Tb = T;                              % initial and boundary temperatures
alpha = 0;                                  % heat transfer coefficient
rmax1 = 0; rmax2 = 0; r0 = 0;               % explosion criteria
k = 1; dt0 = 2e-4;                          % time grid step
t = 0; t_end = 1e-1; aI1 = 0; aI2 = 0;      % time and explosion checkmarks
t_ign1 = 0; T_ign1 = 0;                     % variables to store ignition point
while t(end) < t_end                        % time starts
    k = k + 1;                              
    Ts = T(k-1); Ns = N(k-1,:);             % I use semi-implicit methods
    K0 = exp(Ts/(1+Ar*Ts));                 % chemical reaction rate is determined for the previous temperatures
    Ks = exp(-s/Ar).*exp(s*Ts./(1+Ar.*Ts)); % the factors for lower and higher activation energies
    dt = dt0;                               % there should be a way to access is the time step reasonable
    N(k,:) = Ns./(1 + dt*K0*Ks);            % reaction of pseudocomponents
    NS(k) = sum(N(k,:));                    % overall convesion
    if k == 2                               % 
        r0 = (NS(k)-NS(k-1))/dt;            % average conversion rate
    end
    if isothermal                           
        T(k) = T(k-1);                      % temperature does not change
    else
        T(k) = (T(k-1) + (Tm-s)*(N(k-1,:)-N(k,:))')/(1 + dt*alpha); % temperature changes
    end
    smax(k) = s(N(k,:) == max(N(k,:)));     % position of PDF maximum
    t(k) = t(k-1) + dt;                     % the end of time step 
    if aI1 == 0                                             % if explosion has not happen
        if k > 2                                            % then we look at the heating rate
            if (T(k)-T(k-1))/(t(k)-t(k-1)) > rmax1          % and check is it high
                t_ign1 = t(k); T_ign1 = T(k); aI1 = 0;      % record current maximum point
                rmax1 = (T(k)-T(k-1))/(t(k)-t(k-1));        % refresh maximum value
            end
        end
    end
    if aI2 == 0                                             % the same 
        if k > 2                                            % but for the maximum 
            if (NS(k)-NS(k-1))/(t(k)-t(k-1)) > rmax2        % conversion rate
                t_ign2 = t(k); T_ign2 = Tb(k); aI2 = 0;     % record
                rmax2 = (NS(k)-NS(k-1))/(t(k)-t(k-1));      % refresh
            end
        end
    end
    if (NS(k)) < 1e-3;                                      % if conversion is full enough
        t_end = t(end);                                     % then stop copmutations before t_end
    end
    g(k,:) = N(k,:)/NS(k);                                  % current distribution of pseudocomponents
end

Tmax = max(T);                                              % overall temperature maximum
% Plotting graphs
figure
subplot(1,2,1)
plot(t,NS)                                                  % overall conversion
title('Overall conversion')                                 % 
xlabel('Time')                                              %    
ylabel('N/N_0')                                             % 
subplot(1,2,2)
plot(t,T,t_ign1,T_ign1,'o')                                 % temperature dynamics (with ignition point)
title('Temperature')                                        % 
xlabel('Time')
ylabel('\theta')

figure
subplot(1,2,1)
plot(s,N)                                                   % pseudocomponents distribution evolution
axis([-3*sD, 3*sD, 0, max(max(N))])                         % look at 3 sigma
title('Distibution evolution')
xlabel('Deviation')
ylabel('Number')
subplot(1,2,2)
plot(NS, T)                                                 % phase trajectory in conversion-temperature space
title('Phase plane')
ylabel('Temperature')
xlabel('Concentration')

M = [Tm, sD, alpha, Ar, t_ign1, T_ign1, Tmax, NS(end), smax(end) ];     % characterization of the solution
