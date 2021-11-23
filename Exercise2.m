%=============================================
% Exercise 2: in this exercise, we implement the Dirichlet-Neumann method for
% the solution of a PDE described in the file room_data.
%=============================================
clear all; close all;

room_data;                                                                 % include problem parameters
a=8;                                                                       % location of the interface.
theta=0.5;                                                                 %relaxation parameter
pe=1e12;                                                                   % large Robin parameter to emulate a Dirichlet condition by penalty                                                                         
A=A2d(eta,h,J+2,J);                                                        % global stiffness matrix
u=Solve2dR(A,f,h,J+2,J,pe*gg,pe*gd,pe,pe);                                 % Solve global problem

f1=f(:,2:a);                                                               % in \Omega_1 we solve a complete Dirichlet Problem.
f2=f(:,a+1:end);                                                           % in \Omega_2 we solve a Neumann(left) Dirichlet(right- using penalization) problem.
g=zeros(J,1);
x1=(0:h:a*h);
x2=(a*h:h:1);
y=(0:h:1);
z1=zeros(1,a+1);z2=zeros(1,J-a+2);                                         % for plotting purposes
e=ones(J,1);
Na=[sparse(eye(J,J)),-sparse(diag(-e(1:end-1)/2,-1)+diag((eta*h^2+4)*e/2)+diag(-e(1:end-1)/2,1))]/h; %operator to extract Neumann data.
errorvet(1)=norm(u,2);
Nx1=a-1;  Nx2=J+2-a;                                                       % Number of points in each subdomain
Ny1=J; Ny2=J;
A1=A2d(eta,h,Nx1,Ny1);                                                     % subdomain matrices
A2=A2d(eta,h,Nx2,Ny2);

for i=1:maxiter                                                            % start iterations
    u1=Solve2d(A1,f1,h,Nx1,Ny1,gg,g);                                      % Solve Dirichlet problem
    ta=Na*[u1(:,end-1);u1(:,end)]+f2(:,1)*h/2;                             % compute Neumann derivative 
    u2=Solve2dR(A2,f2,h,Nx2,Ny2,ta,gd,0,pe);                               % Solve Neumann problem
    g=theta*g+(1-theta)*u2(:,1);                                           % update the trace
    ufin=[u1(:,1:a),(u1(:,a+1)+u2(:,1))/2,u2(:,2:end)];
    errorvet(i+1)=norm(u-ufin,2);                                          % compute error
    if plt==1
        mesh(x1,y,[z1;u1;z1]); hold on; mesh(x2,y,[z2;u2;z2]);hold off;
        xlabel('x');ylabel('y');zlabel('Dirichlet-Neumann iterates');
        pause %matlab
        %pause(0.5) %octave
end
end


%=== Verify convergence factor
% figure(2)
% semilogy(iters,errorvet)
% hold on;
% alpha=(a)*h;
% beta=(J-a+1)*h;
% k=(1:J)*pi;
% rho=max(abs(theta-(1-theta)*(tanh(k*beta)./tanh(k*alpha))));
% semilogy(iters,abs(rho).^(iters))
% grid on
% xlabel('Iterations');
% ylabel('Error');
% legend('Error','Theoretical decay')
