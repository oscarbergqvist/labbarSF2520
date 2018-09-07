function u_t=f5(u)
k=0.065;
u_t(1)= -k*u(1)* ((u(1)^2+u(3)^2)^(1/2));
u_t(2)= u(1);
u_t(3)= -9.82 - k*u(3)*(u(1)^2+u(3)^2)^(1/2); 
u_t(4)= u(3);
u_t=transpose(u_t);  
end