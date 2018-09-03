function u = f(u1, u2)
    u1_t = (1-u2^2) * u1 - u2;
    u2_t = u1;
    u = [u1_t, u2_t];
end