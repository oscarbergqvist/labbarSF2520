%% PART 1

% BC = "sin";
BC = "square";

<<<<<<< HEAD
BC = 2; % sätts till nummer av rand värde
x_max= 2; tau=1; t_max= 2*tau; N=100;  a=1;
sigm = 0.8; % sigma is the corrant number
% startar på 0 i index och går till
=======
sigma = 0.9;     
>>>>>>> bacc90d71d69a2e103cfd6771a07fc49012c0cc4

% sigma is the courant number
% All three methods are stable for 0 < sigma <= 1 
% This is in accorance with the theoretical results 
% We have the same stability region for square and sine BC 
% The Upwind and Lax Fredrich methods gives smeared numerical solutions 
% The The Lax Wendroff introduce oscillations for square BC 
% The best method for the sin BC is the upwind method 
% The best method for the square BC is the lax wendroff

[u0, x_vec, t_vec] = get_u0(BC, sigma);

N = length(x_vec) - 1;
M = length(t_vec) - 1;

% Upwind
u_upwind = upwind(u0, sigma, N, M);

% Lax-Fredrich
u_fred = lax_fredrich(u0, sigma, N, M);

% Lax Wendroff
u_wen = lax_wendroff(u0, sigma, N, M);

<<<<<<< HEAD
sigm = 0.1;
h_x = 0.1;
=======
>>>>>>> bacc90d71d69a2e103cfd6771a07fc49012c0cc4

% figure('Name', 'Upwind')
% plot(x_vec, u_upwind(:, end));

% figure('Name', 'LaxFredrich')
% plot(x_vec, u_fred(:, end));

% figure('Name', 'LaxWendroff')
% plot(x_vec, u_wen(:, end));


figure('Name', 'Upwind, LaxFredrich, LaxWendroff')
up_plot = plot(x_vec, u_upwind(:, end));
hold on
fred_plot = plot(x_vec, u_fred(:, end));
hold on
wen_plot=plot(x_vec, u_wen(:, end));
legend([up_plot, fred_plot, wen_plot], ...
        "Upwind", "Lax Freidrich", "Lax Wendroff")


%% Part 2a 

sigma = 1.0; 
dx = 0.05;

[u0, x_vec, t_vec] = get_u0_part2(sigma, dx);
N = length(x_vec) - 1;
M = length(t_vec) - 1;

u_upwind = upwind_part2(u0, sigma, N, M);
u_wen = lax_wendroff_part2(u0, sigma, N, M);

figure('Name', 'Upwind part 2')
subplot(1,2,1)
surf(t_vec, x_vec, u_upwind)
subplot(1,2,2)
surf(t_vec, x_vec, u_wen)

%% PART 2b

% The upwind method smooths out the solutions 
% The more accurate method is the lax wendroff method 
% Higher sigma gives less smoothing for upwind
%
% sigma = 1.0;
% dx = 0.01;
% sigma = 0.1;
% dx = 0.01;
% sigma = 0.8;
% dx = 0.01;
sigma = 0.1;
dx = 0.01;
% sigma = 0.5;
% dx = 0.2;

[u0, x_vec, t_vec] = get_u0_part2(sigma, dx);
N = length(x_vec) - 1;
M = length(t_vec) - 1;

u_upwind = upwind_part2(u0, sigma, N, M);
u_wen = lax_wendroff_part2(u0, sigma, N, M);

figure('Name', 't = 2.5')
up_plot_25 = plot(x_vec, u_upwind(:, floor((M+1)/2)));
hold on 
wen_plot_25 = plot(x_vec, u_wen(:, floor((M+1)/2)));
legend([up_plot_25, wen_plot_25], ...
    'Upwind at t = 2.5', 'Wendroff at t = 2.5')
hold off

figure('Name', 't = 5')
up_plot_5 = plot(x_vec, u_upwind(:, floor(M+1)));
hold on 
wen_plot_5 = plot(x_vec, u_wen(:, floor(M+1)));
legend([up_plot_5, wen_plot_5], 'Upwind at t = 5', 'Wendroff at t = 5')


%%

function [u, x_vec, t_vec] = get_u0(BC, sigm)
 
    tau=1; 
    a=1;  
    x_max= 2; 
    t_max= 2*tau;     
    N=100;

    h_x = x_max/N;
    h_t = (sigm/a) *  h_x;  
    t_vec = 0:h_t:t_max;
    x_vec = 0:h_x:x_max;    
    
    u = zeros(N+1,length(t_vec));
    
    % boundary conditions
    if BC == "sin"
        u(1,:) =  sin(2*pi*t_vec/tau);
    elseif BC == "square"
        for t = 1:(t_max/h_t+1)
            if mod(ceil(2*(t-1)*h_t/tau), 2) == 0  
                u(1, t) = -1;
            else 
                u(1, t) = 1;
            end
        end
    else 
        error("wrong bc, select 'sin' or 'square'")
        return
    end
end


function u = upwind(u0, sigm, N, M)
    u = u0;
    for t_i = 1:M
        for x_i = 2:(N+1)
            u(x_i,t_i+1)=(1-sigm)*u(x_i,t_i)...
                + sigm*u(x_i-1,t_i);
        end
    end
end
<<<<<<< HEAD
figure
surf(t_vec,x_vec,u_upwind)

=======

function u = lax_fredrich(u0, sigm, N, M)
    u = u0;
    for t_i = 1:M
        for x_i = 2:N
            u(x_i,t_i+1)=(u(x_i-1,t_i)+u(x_i+1,t_i))/2 ...
                -(sigm/2)*(u(x_i+1,t_i) -u(x_i-1,t_i));
        end
        u(N+1,t_i+1)= 2*u(N,t_i+1)-u(N-1,t_i+1); 
    end
end
>>>>>>> bacc90d71d69a2e103cfd6771a07fc49012c0cc4

function u = lax_wendroff(u0, sigm, N, M)
    u = u0;
    for t_i = 1:M
        for x_i = 2:N
            u(x_i,t_i+1)= u(x_i,t_i) ...
                -(sigm/2)*(u(x_i+1,t_i) -u(x_i-1,t_i))...
                +sigm^2/2*(u(x_i+1,t_i) -2*u(x_i,t_i)+u(x_i-1,t_i));
        end
        u(N+1,t_i+1)= 2*u(N,t_i+1)-u(N-1,t_i+1); 
    end
end

function [u, x_vec, t_vec] = get_u0_part2(sigm, h_x) 

    T_cool = 10; 
    T_hot = 100; 
    L = 3; 
    k = 0.2; 
    v = 1;
    t_max = 5; 
    x_max = L;
    
    N = x_max / h_x;
    h_t = (sigm/v) *  h_x;
    M = t_max / h_t;

    t_vec = 0:h_t:t_max;
    x_vec = 0:h_x:x_max;

    u = zeros(N+1, M+1);
    u(:,1)= T_cool;

    for i=1:M+1
        if (i-1) * h_t <= 0.5
            u(1,i) = T_cool + (T_hot-T_cool) * sin(pi*(i-1)*h_t);
        elseif (i-1) * h_t<=2
            u(1,i) = T_hot;
        else
            u(1,i) = T_hot + T_cool * sin(4*pi*((i-1)*h_t-2));    
        end
    end
end
<<<<<<< HEAD
figure
surf(t_vec,x_vec,u_wen)

% %b)
figure(1)
up_plot_25 = plot(x_vec, u_upwind(:, floor(2.5/h_t + 1)));
hold on 
wen_plot_25 = plot(x_vec, u_wen(:, floor(2.5/h_t + 1)));
legend([up_plot_25, wen_plot_25], 'Upwind at t = 2.5', 'Wendroff at t = 2.5')
hold off
=======

function u = upwind_part2(u0, sigm, N, M)
    k = 0.2; 
    t_max = 5; 
    T_cool = 10; 
    u = u0;
    for t_i = 1:M
        for x_i = 2:(N+1)
            u(x_i,t_i+1) = (1-sigm) * u(x_i,t_i) ...
                + sigm*u(x_i-1,t_i) - t_max/M * (k*u(x_i,t_i) - k*T_cool);
        end
    end
end
>>>>>>> bacc90d71d69a2e103cfd6771a07fc49012c0cc4

function u = lax_wendroff_part2(u0, sigm, N, M)
    k = 0.2; 
    t_max = 5; 
    T_cool = 10; 
    h_t = t_max/M;
    u = u0;
    for t_i = 1:M
        for x_i = 2:N
            u(x_i,t_i+1) = u(x_i,t_i) - sigm * (1-k*h_t)/2 ...
                * (u(x_i+1,t_i) - u(x_i-1,t_i)) ...
                + sigm^2/2*(u(x_i+1,t_i)-2*u(x_i,t_i)+u(x_i-1,t_i))...
                - h_t * (1-k*h_t/2) * (k*u(x_i,t_i)-k*T_cool);
        end
        u(N+1,t_i+1) = 2*u(N,t_i+1)-u(N-1,t_i+1); 
    end
end