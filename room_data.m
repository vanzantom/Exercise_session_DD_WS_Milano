% Parameters of the problem
eta=0; %reaction term
J=20;  %number of interior mesh points in x and y direction
h=1/(J+1); %mesh size
x=(0:1/(J+1):1); y=x; %finite difference mesh, including boundary
f=zeros(J,J+2);  % source term include boundary
xi=x(2:end-1); yi=xi; %finite difference mesh, excluding boundary
f([yi>0.4 & yi<0.6],[xi>0.4 & xi<0.6])=50; %define force term
gg=0.3*ones(J,1); gg(yi>0.5 & yi<0.9)=1; %boundary data.
gd=zeros(J,1);

%Others parameters:
maxiter=10;%maximum number of iterations
iters=(1:maxiter+1)';
plt=1;% plt=1 to show the subdomain solutions.