function u=Solve2dR(A,f,h,nx,ny,gg,gd,p1,p2)
% SOLVE2dR solves a boundary value problem represented by the matrix A and
% force term f.
% The domain is discretized with Nx\times Ny points, with mesh size h.
% On the horizontal edges, homogeneous boundary conditions are imposed.
% On the left vertical edge, we impose the Robin boundary condition
% \partial_n u+ p1 u =gg, while on the right one,
% \partial_n u+p2 u=gd.
   A(1:ny,1:ny)=A(1:ny,1:ny)/2+p1/h*speye(ny,ny); %
   A(end-ny+1:end,end-ny+1:end)=A(end-ny+1:end,end-ny+1:end)/2+p2/h*speye(ny,ny);
   f(1:ny,1)=f(1:ny,1)/2+gg/h;                                                  % add boundary conditions into rhs
   f(1:ny,end)=f(1:ny,end)/2+gd/h;
   u=A\f(:);
   u=reshape(u,ny,nx);
end