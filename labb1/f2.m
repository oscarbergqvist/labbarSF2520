function x_t = f2(x)
    k1 = 0.04;
    k2 = 10^4;
    k3 = 3*10^7;
    x_t = zeros(3,1);
    x_t(1) = -k1 * x(1) + k2 * x(2) * x(3);
    x_t(2) = k1 * x(1) - k2 * x(2) * x(3) - k3 * x(2)^2;
    x_t(3) = k3 * x(2)^2;
end 