function U = rungeKutta4(u_0,N,t_max,f_h)

    h = t_max/N;
    t0 = 0;
    t = zeros(N, 1);
    U(:, 1) = vertcat(u_0, t0);
    u=u_0; 
    for k = 1:N
        k1 = f_h(u);
        k2 = f_h(u + h * k1);
        k3 = f_h(u + h*k1/4 + h*k2/4);
        u = u + h/6 * (k1 + k2 + 4*k3);
        t(k+1) = t(k) + h; 
        U(:, k+1) = vertcat(u, t(k+1));
    end
end 