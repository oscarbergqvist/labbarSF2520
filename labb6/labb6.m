% 
% % controllerar att kraven är uppfylda sym och pos def
% % A var lika med sitt conjugat se A==A' nedan 
% isequal(A,A')
% % Är A positiv definit? Testa med commandot chol
% [~,p]=chol(A) % P är noll om A positivt definit


%% set_up
S1=load('cooling_flange.mat');
A_varm = S1.A;
S2=load('convdiff.mat');
A_confdiff= S2.A;

%% Jacobian and CGD
start("CGD",100,[20,1;10,3]')

%% b)
n_vec=8:15;
begin= n_vec(1)-1;
d_vec=1:4;
diff_mat= zeros(length(n_vec),length(d_vec));
for d=d_vec
    for n=n_vec 
        
        N=n^d
        A = lap(n,d);
        b=rand(N,1);    
        x= ones(N,1);
        r_old = b-A*x;
        p= r_old;
        loops=0;
        tic
        while norm(A*x-b)/norm(b)>1e-10
          alpha= r_old'*r_old/(p'*A*p);
          x=x+alpha*p;
          r_new= r_old - alpha*A*p;
          norm(r_new);
          if r_new<1e-10
              break
          end
          beta= r_new'*r_new/(r_old'*r_old);
          p = r_new + beta*p;
          loops= loops +1; 
          r_old=r_new; 

        end
        egen_tid=toc

        tic
        ans=A\b;
        mat_tid=toc
        diff_norm=norm(ans - x);
        diff_mat(n-begin,d)=egen_tid-mat_tid;
    end
end
figure
surf(d_vec,n_vec,diff_mat)
xlabel('d')
ylabel('n')
%% part 2 a)
spy(A_varm)
B_varm= rand(size(A_varm,1),1);
tic
[X,FLAG,RELRES,ITER,RESVEC]=pcg(A_varm,B_varm,1e-4,1e5);
pcg_time= toc
tic
back=A_varm\B_varm;
back_time=toc 
X_diff=norm(back-X)/length(X)
semilogy(RESVEC)

%% 2 b) 
M= diag(diag(A_varm));
B_varm= rand(size(A_varm,1),1);
tic
[X,FLAG,RELRES,ITER,RESVEC]=pcg(A_varm,B_varm,1e-4,1e5,M);
pcg_time= toc
tic
back=A_varm\B_varm;
back_time=toc 
X_diff=norm(back-X)/length(X)
plot(RESVEC)
%% 2 b) 

L = ichol(A_varm);
B_varm= rand(size(A_varm,1),1);
tic
[X,FLAG,RELRES,ITER,RESVEC]=pcg(A_varm,B_varm,1e-4,1e5,L,L');
pcg_time= toc
tic
back=A_varm\B_varm;
back_time=toc 
X_diff=norm(back-X)/length(X)
plot(RESVEC)

%% 2 c) 0
M= diag(diag(A_confdiff));
B_varm= rand(size(A_confdiff,1),1);
tic
[X,FLAG,RELRES,ITER,RESVEC]=pcg(A_confdiff,B_varm,1e-4,1e5,M);
pcg_time= toc
tic
back=A_confdiff\B_varm;
back_time=toc 
X_diff=norm(back-X)/length(X)
plot(RESVEC)
%% 2 c) 1
M= diag(diag(A_confdiff));
B_varm= rand(size(A_confdiff,1),1);
tic
[X,FLAG,RELRES,ITER,RESVEC]=gmres(A_confdiff,B_varm, [], 1e-4,1e5,M);
pcg_time= toc
tic
back=A_confdiff\B_varm;
back_time=toc 
X_diff=norm(back-X)/length(X)
plot(RESVEC)
%% 2 c) 2

[L, U] = ilu(A_confdiff);
B_varm= rand(size(A_confdiff,1),1);
tic
[X,FLAG,RELRES,ITER,RESVEC]=gmres(A_confdiff,B_varm, [],1e-4,1e5,L,U);
pcg_time= toc
tic
back=A_confdiff\B_varm;
back_time=toc 
X_diff=norm(back-X)/length(X)
plot(RESVEC)
%% Functions 

function start(method,iter_max,nd_seq)
figure 

    for pair=nd_seq    
        pair
        if method == "Jacobian"
            rel_res_seq=jacobian(iter_max,pair);
        elseif method == "CGD"
            rel_res_seq=cgd(iter_max,pair);
        else
            "fel sträng"
        end
        plot(rel_res_seq)
        hold on
%         legend("n,d="+string(pair(1))+","+string(pair(2)))
    end
    legenden = string(nd_seq)';
    legend("n,d="+legenden(:,1)+", "+legenden(:,2));
end

function  rel_res_seq= jacobian(iter_max,pair)

    n=pair(1);
    d=pair(2);
    N=n^d;

    rel_res_seq= zeros(iter_max,1);
    x=zeros(N,1);
    A = lap(n,d);

    b=rand(N,1);
    for i=1:N
        D(i,i)= A(i,i);
    end
    % Jacobi method 
    norm(full(D^(-1)*(A-D)))

    for i=(1:iter_max)
        x=D^(-1)*(b-(A-D)*x);
        rel_res_seq(i)=norm(A*x-b)/norm(b);
    end
end

function rel_res_seq= cgd(iter_max,pair)

	n=pair(1);
    d=pair(2);
    N=n^d;
    A = lap(n,d);
    b=rand(N,1);    
    rel_res_seq= zeros(iter_max,1);
    x= ones(N,1);
    r_old = b-A*x;
    p= r_old;
    loops=0;
    for i=1:iter_max
      alpha= r_old'*r_old/(p'*A*p);
      x=x+alpha*p;
      r_new= r_old - alpha*A*p;
      norm(r_new)
      if r_new<1e-10
          break
      end
      beta= r_new'*r_new/(r_old'*r_old);
      p = r_new + beta*p;
      loops= loops +1; 
      r_old=r_new; 
      rel_res_seq(i)=norm(A*x-b)/norm(b);
    end
end