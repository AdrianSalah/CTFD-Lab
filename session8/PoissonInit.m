function [A] = PoissonInit(N, dx)
%POISSONINIT returns Coefficient matrix for 1D Poisson equation
%   Parameters
%   N:  number of points (of subdomain)
%   dx: spacial resolution

A = zeros(N);

for i=2:N-1
   A(i,i-1) = -1/(dx*dx);
   A(i,i) = 2/(dx*dx);
   A(i,i+1) = -1/(dx*dx);  
end

A(1,1) = 1;
A(N,N) = 1;


end

