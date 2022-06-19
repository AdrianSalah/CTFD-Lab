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

h1 = 4;
hm = 4;            % only necessary for quatratic option 
h2 = 2;
l = 4;

% Number of degrees of freedom (number of nodes per length)

dimX = 35;
dimY = 35;


% Parameter for Conjugated Heat Transfer (For Session 04)
alpha = 2;
Tinf = 90;

% Boundary conditions (Only Dirichlet applied in Session 03) 

boundary.south = 'Dirichlet';
boundary.north = 'Dirichlet';
boundary.east = 'Dirichlet';
boundary.west = 'Dirichlet';

% Values for Dirichlet BC

TD.north = 10;
TD.south=50;
TD.west= 50;
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



