%PART 1

%definiendum & definiens

BC = 2; % sätts till nummer av rand värde
x_max= 2; tau=1; t_max= 2*tau; N=100;  a=1;
sigm = 0.8; % sigma is the corrant number
% startar på 0 i index och går till

h_x = x_max/N;
h_t = (sigm/a) *  h_x;  

t_vec = 0:h_t:t_max;
x_vec = 0:h_x:x_max;

u = zeros(N+1,length(t_vec));

% boundary conditions

if BC == 1
    u(1,:) =  sin(2*pi*t_vec/tau);
elseif BC == 2
    for t = 1:(t_max/h_t+1)
        if mod(ceil(2*(t-1)*h_t/tau), 2) == 0  
            u(1, t) = -1;
        else 
            u(1, t) = 1;
        end
    end
else 
    "fel initail krav valt"
    return
end

% upwind
u_upwind= u ;

for t_i = 1:(t_max/h_t)
    for x_i = 2:(N+1)
        u_upwind(x_i,t_i+1)=(1-sigm)*u_upwind(x_i,t_i)...
            + sigm*u_upwind(x_i-1,t_i);
    end
end


figure(1)
up_plot=plot(x_vec,u_upwind(:,end));
%plot(x_vec,u(1,:))
hold on



% Lax-Fredrich
%första index är rum andra är tid
u_fred = u;

for t_i = 1:(t_max/h_t)
    for x_i = 2:N
        u_fred(x_i,t_i+1)=(u_fred(x_i-1,t_i)+u_fred(x_i+1,t_i))/2 ...
            -(sigm/2)*(u_fred(x_i+1,t_i) -u_fred(x_i-1,t_i));
    end
    u_fred(N+1,t_i+1)= 2*u_fred(N,t_i+1)-u_fred(N-1,t_i+1); 
end

fred_plot=plot(x_vec, u_fred(1:(N+1),end));
hold on

% Lax Wendroff
u_wen=u;

for t_i = 1:(t_max/h_t)
    for x_i = 2:N
        u_wen(x_i,t_i+1)= u_wen(x_i,t_i) ...
            -(sigm/2)*(u_wen(x_i+1,t_i) -u_wen(x_i-1,t_i))...
            +sigm^2/2*(u_wen(x_i+1,t_i) -2*u_wen(x_i,t_i)+u_wen(x_i-1,t_i));
    end
    u_wen(N+1,t_i+1)= 2*u_wen(N,t_i+1)-u_wen(N-1,t_i+1); 
end

wen_plot=plot(x_vec, u_wen(1:(N+1),end));
legend([up_plot, fred_plot, wen_plot],"upwind","Freidrich","Wendroff")


%% Part 2 !!! dododo !!!

% definiens definiedum
T_cool=10; T_hot=100; L=3; k=0.2; v=1;
t_max=5; x_max=L;

% The optimal solution 
% sigm = 1.0;
% h_x = 0.01;

sigm = 0.1;
h_x = 0.1;

% sigm = 0.8;
% h_x = 0.01;

% sigm = 0.8;
% h_x = 0.1/8;

% sigm = 0.2;
% h_x = 0.01;

% sigm = 0.5;
% h_x = 0.01;


N = x_max / h_x;
h_t = (sigm/v) *  h_x;
M = t_max / h_t;

t_vec = 0:h_t:t_max;
x_vec = 0:h_x:x_max;

u = zeros(N+1,M+1);
u(:,1)= T_cool;
% a)

for i=1:M+1
    if (i-1)*h_t<=0.5
        u(1,i)=T_cool+(T_hot-T_cool)*sin(pi*(i-1)*h_t);
    elseif (i-1)*h_t<=2
        u(1,i)=T_hot;
    else
        u(1,i)=T_hot+ T_cool*sin(4*pi*((i-1)*h_t-2));    
    end
end
% Upwind

u_upwind= u ;

for t_i = 1:(t_max/h_t)
    for x_i = 2:(N+1)
        u_upwind(x_i,t_i+1)=(1-sigm)*u_upwind(x_i,t_i)...
            + sigm*u_upwind(x_i-1,t_i)-h_t*(k*u_upwind(x_i,t_i)...
            -k*T_cool);
    end
end
figure
surf(t_vec,x_vec,u_upwind)


%Wendroff
u_wen=u;

for t_i = 1:(t_max/h_t)
    for x_i = 2:N
        u_wen(x_i,t_i+1)=u_wen(x_i,t_i) - sigm*(1-k*h_t)/2*(u_wen(x_i+1,t_i)...
            -u_wen(x_i-1,t_i))+ sigm^2/2*(u_wen(x_i+1,t_i)...
            -2*u_wen(x_i,t_i)+u_wen(x_i-1,t_i))...
            -h_t*(1-k*h_t/2)*(k*u_wen(x_i,t_i)-k*T_cool);
    end
    u_wen(N+1,t_i+1)= 2*u_wen(N,t_i+1)-u_wen(N-1,t_i+1); 
end
figure
surf(t_vec,x_vec,u_wen)

% %b)
figure(1)
up_plot_25 = plot(x_vec, u_upwind(:, floor(2.5/h_t + 1)));
hold on 
wen_plot_25 = plot(x_vec, u_wen(:, floor(2.5/h_t + 1)));
legend([up_plot_25, wen_plot_25], 'Upwind at t = 2.5', 'Wendroff at t = 2.5')
hold off

figure(2)
up_plot_5 = plot(x_vec, u_upwind(:, floor(5/h_t + 1)));
hold on 
wen_plot_5 = plot(x_vec, u_wen(:, floor(5/h_t + 1)));
legend([up_plot_5, wen_plot_5], 'Upwind at t = 5', 'Wendroff at t = 5')
