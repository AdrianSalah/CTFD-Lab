function [b] = PoissonRHS(N, dx, xs)
%POISSONRHS returns RHS vector for 1D poisson problem
%   Parameters:
%   N:      number of points (of subdomain)
%   dx:     spacial resolution
%   xs:     starting point

b = zeros(N,1);
xx = xs:dx:(N-1)*dx + xs;
b = sin(pi*xx); 
b(1) = 0; b(end) = 0;
end

