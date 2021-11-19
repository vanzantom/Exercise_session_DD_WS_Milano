function [A]=A2d(eta,h,Nx,Ny)
% A2D finite difference approximation of eta-Delta in 2d
% A =A2d(eta,Nx,Ny,h); constructs the finite difference approximation
% to eta-Delta on a Nx x Ny grid with spacing h=1/(Nx+1) both in x and y, and homogeneous
% Dirichlet conditions all around
% Construct 1D Laplace operator in x
Ax=(1/h^2)*spdiags([-ones(Nx,1) (2)*ones(Nx,1) -ones(Nx,1)],-1:1,Nx,Nx);
% Construct 1D Laplace operator in y
Ay=(1/h^2)*spdiags([-ones(Ny,1) 2*ones(Ny,1) -ones(Ny,1)],-1:1,Ny,Ny);
% Construct the 2D Laplace operator in 2D (with Dirichlet b.c. and per subdomain)
A=kron(Ax,speye(Ny))+kron(speye(Nx),Ay);
A=eta*speye(Ny*Nx)+A;

end