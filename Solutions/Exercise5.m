%=============================================
% Exercise 5: in this exercise, you implement the Neumann-Neumann method for
% the solution of a PDE described in the file room_data.
%=============================================
clear all; close all;

room_data;                                     
pe=1e12;                                                                   % penalization parameter to impose weakly Dirichlet B.C.
a=6;                                                                       % location of interface
th=0.25;                                                                    %relaxation parameter
f1=f(:,1:a-1); f2=f(:,a+2:end-1);                                          % local force terms for Dirichlet solves
f_n=zeros(J,J+2);
f1_n=f_n(:,1:a+1); f2_n=f_n(:,a+1:end);                                    % local force terms for Neumann solves
g=zeros(J,1);                                                              % updated trace
gzero=zeros(J,1);                                                          % boundary condition for phi1 and phi2
x1=(0:h:a*h);
x2=(a*h:h:1);
y=(0:h:1);
z1=zeros(1,a+1);z2=zeros(1,J-a+2); z3=zeros(1,J+2);                        % for plotting purposes
e=ones(J,1);
A=A2d(eta,h,J,J);                                                          % global stiffness matrix
u=Solve2d(A,f(:,2:J+1),h,J,J,gg,gd);                                       % Solve global problem
u=[z3;u;z3]; 
Nx1=a-1;  Nx2=J-a;                                                       % Number of points in each subdomain for Dirichlet solve
Ny1=J; Ny2=J;
A1=A2d(eta,h,Nx1,Ny1);                                                     % assemble subdomain stiffness matrices Dirichlet
A2=A2d(eta,h,Nx2,Ny2);

Ny1n=Ny1; Ny2n=Ny2;
Nx1n=a+1;  Nx2n=J+2-a;                                                       % Number of points in each subdomain for Dirichlet solve
A1n=A2d(eta,h,Nx1n,Ny1);                                                     % assemble subdomain stiffness matrices Dirichlet
A2n=A2d(eta,h,Nx2n,Ny2);


errorvet(1)=norm(u,2);
Na=[sparse(eye(J,J)),-sparse(diag(-e(1:end-1)/2,-1)+diag((eta*h^2+4)*e/2)+diag(-e(1:end-1)/2,1))]/h;        %normal derivatives
Nb=[-sparse(diag(-e(1:end-1)/2,-1)+diag((eta*h^2+4)*e/2)+diag(-e(1:end-1)/2,1)),sparse(eye(J,J))]/h;
for i=1:maxiter
    u1=Solve2d(A1,f1,h,Nx1,Ny1,gg,g);                                           % Dirichlet solve in each subdomain     
    u2=Solve2d(A2,f2,h,Nx2,Ny2,g,gd);
    gn=-Na*[u1(:,end-1);u1(:,end)]-Nb*[u2(:,1);u2(:,2)]-f2(:,1)*h-f1(:,end)*h; %compute difference Neumann data.
    phi1=Solve2dR(A1n,f1_n,h,Nx1n,Ny1n,pe*gzero,gn,pe,0);                          % Neumann solves in each subdomain
    phi2=Solve2dR(A2n,f2_n,h,Nx2n,Ny2n,gn,pe*gzero,0,pe);
    g=g-th*(phi1(:,end)+phi2(:,1));                                        %update trace
    ufin=[u1(:,1:end-1),(u1(:,end)+u2(:,1))/2,u2(:,2:end)];
    errorvet(i+1)=norm(u-[z3;ufin;z3],2);
    if plt==1
        mesh(x1,y,[z1;u1;z1]); hold on;
        mesh(x2,y,[z2;u2;z2]); hold on;
        xlabel('x');ylabel('y');zlabel('Neumann-Neumann iterates');
        pause
    end
end
figure(2)
semilogy((1:maxiter+1)',errorvet)
hold on
alpha=(a)*h;
beta=(J-a+1)*h;
k=(1:J)*pi;
rho=1-th*(tanh(k*alpha)+tanh(k*beta)).*(coth(k*beta)+coth(k*alpha));        % convergence factor seen in the lecture
semilogy(iters,max(abs(rho)).^iters,'r')                                   % plot theoretical convergence rate
grid on
xlabel('Iterations');
ylabel('Error');
legend('Error','Theoretical decay')
