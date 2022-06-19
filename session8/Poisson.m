% 1D Poisson equation

%% set up linear system
L = 1;          % length of domain
N = 13;         % number of points
dx = 1/(N-1);   % spacial resolution
x = 0:dx:1;     % spacial arry

% initialization
A = zeros(N);
b = zeros(1,N);
%sol = zeros(N,1);
% building coefficient matrix

for i=2:N-1
   A(i,i-1) = -1/(dx*dx);
   A(i,i) = 2/(dx*dx);
   A(i,i+1) = -1/(dx*dx);  
end

A(1,1) = 1;
A(N,N) = 1;

b = sin(pi*x);
b(1) = 0;
b(end) = 0;


%% solution using matlab \ operator 
u = A\b';

figure(1)

plot(x, u)
title('Temperature distribution 1D Poisson equation')
xlabel('x')
ylabel('T')


% TODO: how to make color plot of 1D 
% figure(2)

%% assignment 1: parametric 

