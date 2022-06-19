%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the parameters of the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Geometry:
%
%    |---
%    |   ---
%    |      ----
%    |          |
% h1 |----------|  h2  <- symmetry axis
%    |          |
%    |      ----
%    |   ---
%    |---
%
%    |<--  l -->|


% Defin dimension of the trapezoidal domain
% h2 <= h1 !


shape = 'quadratic'  % 'linear' or 'quadratic'

h1 = 10;
hm = 4;            % only necessary for quatratic option 
h2 = 3;
l = 10;

% Number of degrees of freedom (number of nodes per length)

dimX = 10;
dimY = 10;


% Parameter for Conjugated Heat Transfer (For Session 04)
alpha = 5;
betanorth = 0;
betaeast = 4;
betawest = -5;
Tinf = 90;
k = 1;

% Boundary conditions (Only Dirichlet applied in Session 03) 

boundary.south = 'Dirichlet';
boundary.north = 'Dirichlet';
boundary.east = 'Dirichlet';
boundary.west = 'Dirichlet';

% Values for Dirichlet BC

TD.north = 10;
TD.south= 50;
TD.west= 100;
TD.east=10;

%steady or transient behaviour
tempdev = 'transient';
timerange = 1;
tsteps = 80 ;
dt = timerange/tsteps;
time = [timerange, tsteps, dt];

%define preferred solver: 
solver.type = 'theta'; % choose here 'theta' or 'RK4'
solver.type2 = 1; % adjust theta for theta scheme: 0/1 ==> explicit/ implicit, 0.5 Crank-Nicolson, other values [0...1]

% Shape of the Cooling F    in h2 <= h1 !

switch shape
    
    case 'linear'

        formfunction = @(xnorm) (1-xnorm)*h1/2 + xnorm*h2/2;

    case 'quadratic'

        c1 = h2+2*h1/2-2*hm;
        c2 = 2*hm - 3*h1/2 - h2/2;
        c3 = h1/2;

        formfunction = @(xnorm) c1*xnorm.^2 +c2*xnorm + c3;

end



