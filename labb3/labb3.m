%% Forward Euler 

t_p = 1;
N = 10;
K = 400;
T_max = 2;
u0 = zeros(N, 1);
dt = T_max / K;

U = zeros(N+1, K+1);

u = u0;
for i = 1:K+1
    tau = dt*(i-1);
    u_xi0 = double(tau <= t_p);
    U(:, i) = [u_xi0; u];
    u = u + dt * build_u_t(u,tau,N, t_p);
end

surf(linspace(0, 1, K+1), linspace(0, 1, N+1), U)
plot(linspace(0, 1, 11), U(:, 400))
xlabel('$\xi$', 'interpreter', 'latex')
ylabel('u')
axis([0 1 0 1])
%% ODE23
N = 40;

fprintf("ODE23\n")
tic
[TOUT_ode23, YOUT_ode23] = ode23(@f, [0, 2], zeros(N, 1));
toc
timesteps_ode23 = length(TOUT_ode23) - 1
dT_max_ode23 = max(TOUT_ode23(2:end) - TOUT_ode23(1:end-1))


%% ODE23s
N = 40;

A = N^2 * (-2 * diag(ones(N,1), 0) + ...
    diag(ones(N-1, 1), -1) + ...
    diag(ones(N-1, 1), 1)); 
A(N:N-1) = 2;

S = (diag(ones(N,1), 0) + ...
    diag(ones(N-1, 1), -1) + ...
    diag(ones(N-1, 1), 1)); 

fprintf("ODE23s\n")
tic
[TOUT, YOUT] = ode23s(@f, [0, 2], zeros(N, 1), odeset('Jacobian', A));
toc
timesteps = length(TOUT)-1
dT_max = max(TOUT(2:end) - TOUT(1:end-1))


%% FUNCTIONS

function u_t = build_u_t(u, tau, N, t_p) 
    h = 1/N;
    u_t = zeros(N, 1);
    if tau <= t_p
        u_t(1) = 1/h^2 * ((-2 * u(1) + u(2)) + 1); 
    else
        u_t(1) = 1/h^2 * (-2 * u(1) + u(2));
    end
    for i = 2:N-1
        u_t(i) = 1/h^2 * (u(i-1) - 2 * u(i) + u(i+1));
    end
    u_t(N) = 1/h^2 * (2 * u(N-1) - 2 * u(N));
end 


function u_t = f(tau, u) 
    N = 40; 
    h = 1/N;
    u_t = zeros(N, 1);
    if tau <= 1
        u_t(1) = 1/h^2 * ((-2 * u(1) + u(2)) + 1); 
    else
        u_t(1) = 1/h^2 * (-2 * u(1) + u(2));
    end
    for i = 2:N-1
        u_t(i) = 1/h^2 * (u(i-1) - 2 * u(i) + u(i+1));
    end
    u_t(N) = 1/h^2 * (2 * u(N-1) - 2 * u(N));
end 
