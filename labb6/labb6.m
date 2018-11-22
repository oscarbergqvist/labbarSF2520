
%konstanter
n=10;
d=3;

%
N=n^d;
x=zeros(N,1);
A = lap(n,d);

b=rand(N,1);
for i=1:N
    D(i,i)= A(i,i);
end

norm(full(D^(-1)*(A-D)))
x_ny=50*ones(N,1);
x_gammal= x ; 
tol= 10^(-8);
while (norm((x_ny-x_gammal))/norm(x_ny)>tol)
    x_gammal=x_ny;
    
    x_ny=D^(-1)*(b-(A-D)*x_gammal);
end

