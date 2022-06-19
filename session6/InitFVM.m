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


shape = 'quadratic';  % 'linear' or 'quadratic'

h1 = 10;
hm = 4;            % only necessary for quatratic option 
h2 = 3;
l = 10;

% Number of degrees of freedom (number of nodes per length)

dimX = 10;
dimY = 10;


% Parameter for Conjugated Heat Transfer (For Session 04)
alpha = 5;
Tinf = 90;

% specify if steady or unsteady simulation
simulationType = 'steady'
%simulationType = 'transient';

% specify type of time integration
%timeintegrationType = 'explicitEuler'
%timeintegrationType = 'implicitEuler'
timeintegrationType = 'theta'
%timeintegrationType = 'rungeKutta4'

% specify parameter for theta scheme
theta = 1; %  theta = 0 -> explicit euler, theta = 1 -> implicit euler

% Timestepsize and Endtime for unsteady case
dt = 0.1;
tend = 10;

% Boundary conditions (Only Dirichlet applied in Session 03) 

boundary.south = 'Neumann';
boundary.north = 'Robin';
boundary.east = 'Robin';
boundary.west = 'Dirichlet';

% Values for Dirichlet BC

TD.north = 10;
TD.south=50;
TD.west= 100;
TD.east=10;


% Shape of the Cooling Fin

% h2 <= h1 !

switch shape
    
    case 'linear'

        formfunction = @(xnorm) (1-xnorm)*h1/2 + xnorm*h2/2;

    case 'quadratic'

        c1 = h2+2*h1/2-2*hm;
        c2 = 2*hm - 3*h1/2 - h2/2;
        c3 = h1/2;

        formfunction = @(xnorm) c1*xnorm.^2 +c2*xnorm + c3;
        
%     case 'sin'
%         
%         formfunction = @(xnorm) sin(xnorm)./xnorm;

end



