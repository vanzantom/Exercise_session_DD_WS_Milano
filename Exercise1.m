%=============================================
% Exercise 1: in this exercise, we implement an alternating Schwarz method for
% the solution of a PDE described in the file room_dataSchwarz.
%=============================================

clear all;
close all;
room_data;                                                                 % load problem's parameters
a=8; d=4;                                                                  % spatial decomposition: First sub. takes Dirichlet data at x=a+d. Second sub. takes Dirichlet data at x=a                                                        
f=f(:,2:end-1);                                                            % restrict f to the interior of Omega.
f1=f(:,1:a+d-1); f2=f(:,a+1:end);                                          % a+d is boundary term of \Omega_1. a left boundary of \Omega_2
u1=[gg zeros(J,a+d)];                                                      % zero initial guess, and add boundary value
u2=[zeros(J,J-a+1) gd];
y=0:h:1;                                                                   % finite difference meshes
x1=(0:h:(a+d)*h); x2=(a*h:h:1);
Nx1=a+d-1;  Nx2=J-a;                                                       % number of interior points in \Omega_1 and \Omega_2
Ny1=J; Ny2=J;
A1=A2d(eta,h,Nx1,Ny1);                                                     % assemble subdomain stiffness matrices
A2=A2d(eta,h,Nx2,Ny2);
A=A2d(eta,h,J,J);
u=Solve2d(A,f,h,J,J,gg,gd);                                                % solve global problem
z1=zeros(1,a+d+1);z2=zeros(1,J-a+2); z3=zeros(1,J+2);                      % for plotting purposes
u=[z3;u;z3];                                                               % Add horizontal homogeneous Dirichlet Bc.
errorvet=ones(maxiter,1);
errorvet(1)=norm(u,2);


for i=1:maxiter                                                            % start iteration Schwarz iteration
    u1n=Solve2d(A1,f1,h,Nx1,Ny1,gg,u2(:,d+1));                             % solve first subdomain
    u1=u1n;                                                                % update the iterates
    u2n=Solve2d(A2,f2,h,Nx2,Ny2,u1(:,end-d),gd);                           % solve second subdomain
    u2=u2n;
    ufin=[u1n(:,1:a),(u1n(:,a+1:a+d+1)+u2n(:,1:d+1))/2,u2n(:,d+2:end)];    % merge the two contributions in the overlap
    errorvet(i+1)=norm(u-[z3;ufin;z3],2);                                  % compute error between iterate and exact solution                                          
    if plt==1                                                              % plot subdomain solutions
        mesh(x1,y,[z1;u1n;z1]); hold on; mesh(x2,y,[z2;u2n;z2]);hold off;
        xlabel('x');ylabel('y');zlabel('Schwarz iterates');
        pause %matlab
        %pause(0.5) %octave command
    end
end


%=== Verify convergence rate
% figure(2)
% semilogy(iters,errorvet)                                                   % plot error
% hold on
% beta=a*h;                                                                  % location on the boundary of \Omega_2
% alpha=(a+d)*h;                                                             % location on the boundary of \Omega_1
% k=(1:J)*pi;
% rho=max(abs((sinh(k*beta)./sinh(k*(1-beta))).*(sinh(k*(1-alpha))./sinh(k*alpha))));% convergence factor seen in the lecture
% semilogy(iters,rho.^iters,'r')                                             % plot theoretical convergence rate
% grid on
% xlabel('Iterations');
% ylabel('Error');
% legend('Error','Theoretical decay')