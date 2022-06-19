clear all

% Numerical code to solve the 2D steady heat equation by Finite Differences
% Second order accurate in Space
%
%

% Geometrical Parameters

l = 1;                      % Length of the domain in x 
h = 1;                      % Length of the domain in y
dimX = 4;                  % Number of nodes in x
dimY = 4;                  % Number of nodes in y 

% Defining Boundary conditions

% Type
boundary.south = 'Neumann';
boundary.north = 'Robin';
boundary.east  = 'Robin';
boundary.west  = 'Dirichlet';

% Value for Dirichlet BC
TNorth = 10;
TSouth= 50;
TWest= 50;
TEast=10;

% Values of Neumann BC
beta = 0;                   % Heat flux
% Values of Robin BC
alpha = 1;                  % Convective heat transfer coefficient
Tinf = 100;                  % Temperature of the surrounding fluid


% index = @(ii,jj) ii+(jj-1)*dimY;
% Thermal conductivity Coefficient

heat_conduc='homogeneous'   %1) homgeneous, 2) non homogeneous
Kval = 1;


% Region with different K

KnH = 10;                   % Thermal conductivity for the particular region

yk = [0.5 0.8];
xk = [1 1.5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------
% Define Temperature T matrix (Intial values with zero for example)
T = zeros(dimX*dimY, 1);


%------------------------------------------------------
% Create mesh (2D MESH). (Matrix X and Matrix Y)

% for plotting??



%----------------------------------------------------------------------
% Define the Heat conductivity coefficient values


switch heat_conduc
    
    case 'homogeneous'


        
    case 'non homogeneous'


end

%-----------------------------------------------------------------------
% Defining the Source



%-----------------------------------------------------------------------
% Defining Boundary conditions


%------------------------------------------------------------------------
% Constructing Matrix A;
% initialize A


        
%------------------------------------------------------------------------      
% Solving the linear system (use '\' operator)
        

        
%------------------------------------------------------------------------
% Ploting Results

% plot meshgrid (optional)

figure(1)


%surf(x, y, K)

% Do a surface plot
figure(2)


% Do a contour plot

figure(3) 


% Do a color plot





