% labb1
%% PART 1

k_max_vec = [10, 20, 40, 80, 160, 320];

Y = cell(6, 1);
T = cell(6, 1);
i = 1;
for k_max = k_max_vec
    h = 1/ k_max;
    t0 = 0;
    u = [0, 1];
    y = zeros(k_max, 1);
    t = zeros(k_max, 1);
    t(1) = 0;
    y(1) = 1;
    for k = 1:k_max
        k1 = f(u(1), u(2));
        k2 = f(u(1) + h * k1(1), u(2) + h * k1(2));
        k3 = f(u(1) + h*k1(1)/4 + h*k2(1)/4, u(2) + h*k1(2)/4 + h*k2(2)/4);
        u = u + h/6 * (k1 + k2 + 4*k3);
        t(k+1) = t(k) + h; 
        y(k+1) = u(2);
    end
    Y{i} = y
    T{i} = t
    i = i + 1;
    
end 

%%

y = Y{6};
t = T{6};
plot(t, y)

%%

for i = 1:5
    e(i) = Y{i}(end) - Y{6}(end);
    h = 1./k_max_vec;
end

loglog(h(1:5), abs(e))


%% PART 2

k_max_vec = [125, 250, 500, 1000, 2000];

U_cell = cell(5, 1);
i = 1;
for k_max = k_max_vec
    h = 1/ k_max;
    t0 = 0;
    u = [1; 0; 0];
    y = zeros(k_max, 1);
    t = zeros(k_max, 1);
    t(1) = 0;
    U(:, 1) = vertcat(u, t(1));
    for k = 1:k_max
        k1 = f2(u);
        k2 = f2(u + h * k1);
        k3 = f2(u + h*k1/4 + h*k2/4);
        u = u + h/6 * (k1 + k2 + 4*k3);
        t(k+1) = t(k) + h; 
        U(:, k+1) = vertcat(u, t(k+1));
    end
    U_cell{i} = U;
    i = i + 1;
    
end 

%% 

N = 3;
figure('name', 'x1')
loglog(U_cell{N}(4,:), U_cell{N}(1, :))
figure('name', 'x2')
loglog(U_cell{N}(4,:), U_cell{N}(2, :))
figure('name', 'x3')
loglog(U_cell{N}(4,:), U_cell{N}(3, :))


%% ODE23
f2_t_handle = @f2_t;
reltol = [1e-3, 1e-4, 1e-5, 1e-6];
T_OUT = cell(4, 1);
Y_OUT = cell(4, 1);
for i = 1:4
    opts = odeset('RelTol', reltol(i));
    [TO, YO] = ode23(f2_t_handle, [0, 1], [1; 0; 0], opts);
    T_OUT{i} = TO;
    Y_OUT{i} = YO;
    n_steps(i) = length(T_OUT{i}) - 1;
    step_size{i} = T_OUT{i}(2: end) - T_OUT{i}(1:end-1);
end

%%
i = 1
figure('name', '1e-6')
plot(T_OUT{i}(1:end-1), step_size{i})

%% ODE23s
f2_t_handle = @f2_t;
reltol = [1e-3, 1e-4, 1e-5, 1e-6];
T_OUT = cell(4, 1);
Y_OUT = cell(4, 1);
for i = 1:4
    opts = odeset('RelTol', reltol(i));
    [TO, YO] = ode23s(f2_t_handle, [0, 1000], [1; 0; 0], opts);
    T_OUT{i} = TO;
    Y_OUT{i} = YO;
    n_steps(i) = length(T_OUT{i}) - 1;
    step_size{i} = T_OUT{i}(2: end) - T_OUT{i}(1:end-1);
end

%%
i = 4
figure('name', '1e-3')
plot(T_OUT{i}(1:end-1), step_size{i})

