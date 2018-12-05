x_max = 4; 
y_max = 2;
h = 0.1;
u = zeros(N*(M-2), 1);



%% FUNCTIONS 

function [A, b] = gen_system(f, h)
    N = x_max / h - 1;
    M = y_max / h - 1;
    A = sparse(N*M, 1);
    b = zeros(M*N, 1);
    
    % Equations for i = 0
    ind = 1;
    for j = 1:M
        A_i(ind) = (j-1) * N + 1;
        A_j(ind) = (j-1) * N + 1;
        A_v(ind) = 1;
        b((j-1) * N + 1) = 300;
        ind = ind + 1;
    end
    
    % Equations for i = N 
    for j = 1:M
        A_i(ind) = j * N;
        A_j(ind) = j * N;
        A_v(ind) = 1;
        b(j * N) = 600;
        ind = ind + 1;
    end
    
    % Equations for j = 0 
    for i = 2:(N-1)
        A_i(ind) = (i-1);
        A_j(ind) = i;
        A_v(ind) = 1;
        ind = ind + 1;
        A_i(ind) = i;
        A_j(ind) = i;
        A_v(ind) = -4;        
        ind = ind + 1;
        A_i(ind) = i+1;
        A_j(ind) = i;
        A_v(ind) = 1;  
        ind = ind + 1;
        A_i(ind) = i+N;
        A_j(ind) = i;
        A_v(ind) = 2;  
        ind = ind + 1;        
    end
    
    % Equations for j = M 
    for i = (N * (M-1) + 2):(N * M - 1)
        A_i(ind) = (i -1);
        A_j(ind) = i;
        A_v(ind) = 1;
        ind = ind + 1;
        A_i(ind) = i;
        A_j(ind) = i;
        A_v(ind) = -4;        
        ind = ind + 1;
        A_i(ind) = i+1;
        A_j(ind) = i;
        A_v(ind) = 1;  
        ind = ind + 1;
        A_i(ind) = i+N;
        A_j(ind) = i;
        A_v(ind) = 2;  
        ind = ind + 1;     
    end
    
    % Equations for the interior 
    for j = 2:M
        for i = 2:N
            A_i(ind) = (i - 1);
            A_j(ind) = i;
            A_v(ind) = 1;
            ind = ind + 1;
            A_i(ind) = i;
            A_j(ind) = i;
            A_v(ind) = -4;        
            ind = ind + 1;
            A_i(ind) = i+1;
            A_j(ind) = i;
            A_v(ind) = 1;  
            ind = ind + 1;
            A_i(ind) = i+N;
            A_j(ind) = i;
            A_v(ind) = 1;  
            ind = ind + 1;     
        end
    end 
    
end
