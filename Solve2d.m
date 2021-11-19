function u=Solve2d(A,f,h,Nx,Ny,gg,gd)
% SOLVE2D solves a boundary value problem represented by the matrix A and
% force term f.
% The domain is discretized with Nx\times Ny points, with mesh size h.
% On the horizontal edges, homogeneous boundary conditions are imposed.
% On the left vertical edge, we impose u=gg and on the right one, u=gd.

   f(1:Ny,1)=f(1:Ny,1)+gg/h^2;                                             %add boundary conditions into rhs
   f(1:Ny,end)=f(1:Ny,end)+gd/h^2;
   u=A\f(:);
   u=reshape(u,Ny,Nx);
   u=[gg u gd];                                                            %add boundary values to solution
end