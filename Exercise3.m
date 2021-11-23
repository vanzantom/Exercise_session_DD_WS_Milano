%=============================================
% Exercise 3: in this exercise, we implement a parallel optimized Schwarz method for
% the solution of a PDE described in the file room_data.
%=============================================

clear all; close all;
                                                          
room_data;                                                   
a=8; d=4;                                                                  % the right boundary of \Omega_1 is in a+d+1 and the left boundary of \Omega_2 is in a+1
f1=..........;                                                             % Complete!!!
f2=..........;                                               
u1=zeros(J,a+d+1); u2=zeros(J,J-a+2);                                      % zero initial guess
x1=(0:h:(a+d)*h);x2=(a*h:h:1); y=(0:h:1);                                  % finite difference meshes
z1=zeros(1,a+d+1);z2=zeros(1,J-a+2);z3=zeros(1,J+2);                       % for plotting purposes
p=1;                                                                       % nonoptimized parameter
%p=((pi^2+eta)/(d*h))^(1/3);                                               % optimized parameter with overlap.  
%p=sqrt(pi*pi/h);                                                          % optimized parameter without overlap.         
e=ones(J,1);                                                               % support vector to construct normal derivatives
pe=1e12;                                                                   % large Robin parameter to emulate a Dirichlet condition by penalty
A=A2d(eta,h,J+2,J);                                                        % defined global stiffness matrix
u=Solve2dR(A,f,h,J+2,J,pe*gg,pe*gd,pe,pe);                                 % solve global problem
z1=zeros(1,a+d+1);z2=zeros(1,J-a+2); z3=zeros(1,J+2);                      % for plotting purposes
u=[z3;u;z3];                                                               % Add horizontal homogeneous Dirichlet Bc.
Nx1=a+d+1;  
Nx2=.....;                                                                 %Complete!!!
Ny1=J; Ny2=J;
A1=A2d(eta,h,Nx1,Ny1);                                                     % define subdomain matrices
A2=A2d(eta,h,Nx2,Ny2);
errorvet(1)=norm(u,2);                                                  
Na=[sparse(eye(J,J)),-sparse(diag(-e(1:end-1)/2,-1)+diag((eta*h^2+4)*e/2)+diag(-e(1:end-1)/2,1))]/h; %operators which extract Neumann data
Nb=[-sparse(diag(-e(1:end-1)/2,-1)+diag((eta*h^2+4)*e/2)+diag(-e(1:end-1)/2,1)),sparse(eye(J,J)) ]/h;


for i=1:maxiter
    tb=Nb*[u2(:,d+1);u2(:,d+2)]+f2(:,d+1)*h/2+p*u2(:,d+1);                 % derive Robin boundary data from old approximations
    
    ta=................                                                    %Complete!
    
    u1n=Solve2dR(A1,f1,h,Nx1,Ny1,pe*gg,tb,pe,p);                           % Solve subdomain problems
    u2n=..............                                                     %Complete!
    u1=u1n;
    u2=u2n;
    ufin=[u1n(:,1:a),(u1n(:,a+1:a+d+1)+u2n(:,1:d+1))/2,u2n(:,d+2:end)];    % average contribution in the overlap
    errorvet(i+1)=norm(u-[z3;ufin;z3],2);
    if plt==1
        mesh(x1,y,[z1;u1n;z1]); hold on; mesh(x2,y,[z2;u2n;z2]);hold off;  % plot subdomains solves
        pause %matlab
        %pause(0.5) %octave
        xlabel('x');ylabel('y');zlabel('Optimized Schwarz iterates');
   end
end

% %=== Verify convergence factor
% figure(2)                                                                  % plot
% semilogy(iters,errorvet)
% hold on
% 
% %%Complete the code to finish the exercise!
% %beta=......
% %alpha=........
% %k=.........
% %rho=.........
% semilogy(iters,max(rho).^(iters/2),'r')
% grid on
% xlabel('Iterations');
% ylabel('Error');
% legend('Error','Theoretical decay')
