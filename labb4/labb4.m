x_max = 4; 
y_max = 2;
h = 0.1;
u = zeros(N*(M-2), 1);



%% FUNCTIONS 

function [A, b] = gen_system(f, h)
    N = x_max / h - 1;
    M = y_max / h - 1;
    A = sparse(N*(M-2), 1);
    for j = 2:(M+1)
        for i = 2:(N+1)
            if j == 1 || j == M 
                A((M-1) * )
            else
                b(i + j) = -(u(i-1, j) + u(i+1, j) + u(i, j-1) + u(i, j+1) + f);
                A(i+j, i+j) = 
            end
        end
    end
end
