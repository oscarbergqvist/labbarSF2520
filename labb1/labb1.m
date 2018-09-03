% labb1

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


%%



