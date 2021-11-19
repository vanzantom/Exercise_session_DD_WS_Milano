%=============================================
% Exercise 4: in this exercise, you implement the Schwarz method in a substructured form 
% and use GMRES to accederate the convergence of the fixed point iteration.
% We consider again the moded probdem described in room_data
%=============================================
clear all; close all;

room_data;
a=8; d=4;                                                                  % decomposition
f=f(:,2:end-1);                                                            % restrict f to the interior of Omega.
f1=f(:,1:a+d-1);                                                           % local forces
f2=f(:,a+1:end);
x1=0:h:(a+d)*h;                                                            % finite difference meshes
x2=a*h:h:1;
G1=zeros(a+d+1,1);G1(end-d)=1;                                             % construct substructured system:
G2=zeros(J-a+2,1);G2(d+1)=1;
z=zeros(J,1);
Nx1=a+d-1;  Nx2=J-a;                                                       % number of interior points in \Omega_1 and \Omega_2
Ny1=J; Ny2=J;
A1=A2d(eta,h,Nx1,Ny1);                                                     % assemble subdomain stiffness matrices
A2=A2d(eta,h,Nx2,Ny2);

b=[Solve2d(A1,f1,h,Nx1,Ny1,gg,z)*G1;                                       % T*g=b, g unknowns at interfaces
   Solve2d(A2,f2,h,Nx2,Ny2,z,gd)*G2];
T=@(g) [g(1:J)-Solve2d(A1,0*f1,h,Nx1,Ny1,z,g(J+1:2*J))*G1;  
        g(J+1:2*J)-Solve2d(A2,0*f2,h,Nx2,Ny2,g(1:J),z)*G2];

A=A2d(eta,h,J,J);                                                          % compute solution global problem
u=Solve2d(A,f,h,J,J,gg,gd); 
G=zeros(J+2,2);                                                            % assemble operator G which extracts the values at the interior interfaces
G(a+1,1)=1;
G(d+a+1,2)=1;
us=u*G;
us=us(:);

g=zeros(2*J,1);                                                            % zero initial guess
errorvec(1)=norm(b-T(g));                                                  % initial residual
for i=1:maxiter                                                                 % classical Schwarz iteration
  g=g-T(g)+b;                                                              % gn=(I-T)*g+b
  errorvec(i+1)=norm(us-g);                                              % keep residual for plotting
end

semilogy(iters,errorvec)
hold on
beta=a*h;                                                                  % location on the boundary of \Omega_2
alpha=(a+d)*h;                                                             % location on the boundary of \Omega_1
rho=sqrt((sinh(pi*beta)/sinh(pi*(1-beta)))*(sinh(pi*(1-alpha))/sinh(pi*alpha)));% convergence factor seen in the lecture
semilogy(iters,rho.^iters,'r')                                             % plot theoretical convergence rate
xlabel('Iterations');
ylabel('Error/Residual');
grid on
[g,fd,r,it,rk]=gmres(T,b);                % use Krylow solver
semilogy(0:length(rk)-1,rk,'-+');
legend('Error Iterative','Theoretical estimate','GMRES')
