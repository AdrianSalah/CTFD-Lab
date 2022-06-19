clear all

 

% Numerical code to solve the 2D steady heat equation by Finite Differences
% Second order accurate in Space
%
%

 

% Geometrical Parameters
l = 2;                              % Length of the domain in x 
h = 2;                              % Length of the domain in y
dimX = 30;                          % Number of nodes in x
dimY = 30;                          % Number of nodes in y 
index = @(ii,jj) ii+(jj-1)*dimY;    % transforms 2D index into 1D vector index
dimA = dimX*dimY;                   % Dimension of Coefficient matrix A
 

% Defining Boundary conditions
% Type
boundary.south = 'Neumann';
boundary.north = 'Robin';
boundary.east  = 'Robin';
boundary.west  = 'Dirichlet';

% Value for Dirichlet BC
TNorth = 10;
TSouth = 50;
TWest = 50;
TEast = 10;

% Values of Neumann BC
beta = 0;                   % Heat flux
% Values of Robin BC
alpha = 1;                  % Convective heat transfer coefficient
Tinf = 373.15;              % Temperature of the surrounding fluid

% Thermal conductivity Coefficient
heat_conduc = 'non homogeneous'   %1) homgeneous, 2) non homogeneous
Kval = 1;

% Region with different K
KnH = 10;                   % Thermal conductivity for the particular region

yk = [0.5 0.8];             % specifies absolute length of region with different thermal conductivity
xk = [1 1.5];

% Region with source term
s = -1000;
ys = [0.2 0.3];
xs = [0.9 1.1];

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------
% Define Temperature T matrix (Intial values with zero for example)
 
% not used
% Tmat = 273.15 * ones(dimX,dimY);

 

%------------------------------------------------------
% Create mesh (2D MESH). (Matrix X and Matrix Y)

dx = l / (dimX-1);
dy = h / (dimY-1);

 

%----------------------------------------------------------------------
% Define the Heat conductivity coefficient values

% initialize matrix of for thermal conductivity coeff
lam = Kval * ones(dimX,dimY);

switch heat_conduc
    
    case 'homogeneous'
        
    A = zeros(dimA);
    b = zeros(dimA,1);
        
    case 'non homogeneous'
    % assign thermal conductivity for special reagion to thermal conductivity matrix    
    dimKX = floor(xk./l.*dimX);
    dimKY = floor((h-yk)./h.*dimY);

    lam(dimKX(1):dimKX(2),dimKY(2):dimKY(1)) = KnH;

    A = zeros(dimA);
    b = zeros(dimA,1);

end

 

%-----------------------------------------------------------------------
% Defining the Source

smat = zeros(dimX,dimY);

dimsX = floor(xs./l.*dimX);
dimsY = floor((h-ys)./h.*dimY);

smat(dimsY(2):dimsY(1),dimsX(1):dimsX(2)) = s;

b = b + smat(:);
%-----------------------------------------------------------------------
% Defining Boundary conditions

switch boundary.south
    
    case 'Dirichlet'
        i = dimY;
        
        A(dimY:dimY:end,:) = 0;
        
        for j = 1:dimX
            k = index(i,j);
            A(k,k) = 1;
        end

        b(dimY:dimY:end) = TSouth + 273.15;
        
    case 'Neumann'
        i = dimY;
       
        for j = 1:dimX
            k = index(i,j);
            f = -lam(i,j)/(2*dy);
            A(k,k) = A(k,k) - 3*f;
            A(k,k-1) = A(k,k-1) + 4*f;
            A(k,k-2) = A(k,k-2) - 1*f;
        end
        
        
        b(dimY:dimY:end) = b(dimY:dimY:end) + beta;
        
    case 'Robin'
        i = dimY;
        
        for j = 1:dimX
            f = lam(i,j)/(2*dy);
            A(k,k) = A(k,k) + (alpha+3*f);
            A(k,k-1) = A(k,k-1) - 4*f;
            A(k,k-2) = A(k,k-2) + 1*f;   
        end
        
        
        b(dimY:dimY:end) = b(dimY:dimY:end) + alpha*Tinf;
        
end

 

switch boundary.north
    
    case 'Dirichlet'
        i = 1;
        
        A(1:dimY:end,:) = 0;
        
        for j = 1:dimX
            k = index(i,j);
            A(k,k) = 1;
        end
        
        b(1:dimY:end) = TNorth + 273.15;
        
    case 'Neumann'  
        i = 1;
        
        for j = 1:dimX
            k = index(i,j);
            f = -lam(i,j)/(2*dy);
            A(k,k) = A(k,k) - 3*f;
            A(k,k+1) = A(k,k+1) + 4*f;
            A(k,k+2) = A(k,k+2) - 1*f;
            
        end

        b(1:dimY:end) = b(1:dimY:end) + beta;
       
    case 'Robin'    
        i = 1;  
        
        for j = 1:dimX
            k = index(i,j);
            f = lam(i,j)/(2*dy);
            A(k,k) = A(k,k) + (alpha+3*f);
            A(k,k+1) = A(k,k+1) - 4*f;
            A(k,k+2) = A(k,k+2) + 1*f;
        end

        b(1:dimY:end) = b(1:dimY:end) + alpha * Tinf;
        
end

 

switch boundary.east
    
    case 'Dirichlet'      
        j = dimX;
        
        A(end-dimY+1:end,:) = 0;
        
        for i = 1:dimY
            k = index(i,j);
            A(k,k) = 1;
        end
  
        b(end-dimY+1:end) = TEast + 273.15;
        
    case 'Neumann'       
        j = dimX;
 
        for i = 1:dimY
            k = index(i,j);
            f = -lam(i,j)/(2*dx);
            A(k,k) = A(k,k) - 3*f;
            A(k,k-dimY) = A(k,k-dimY) + 4*f;
            A(k,k-2*dimY) = A(k,k-2*dimY) - 1*f;
        end

        
        b(end-dimY+1:end) = b(end-dimY+1:end) + beta;
        
    case 'Robin' 
        j = dimX;
        
        for i = 1:dimY
            k = index(i,j);
            f = lam(i,j)/(2*dx);
            A(k,k) = A(k,k) + (alpha+3*f);
            A(k,k-dimY) = A(k,k-dimY) - 4*f;
            A(k,k-2*dimY) = A(k,k-2*dimY) + 1*f;
        end
            
        b(end-dimY+1:end) = b(end-dimY+1:end) + alpha*Tinf;
end

 

switch boundary.west
    
    case 'Dirichlet'
        j = 1;
        
        A(1:dimY,:) = 0;
        
        for i = 1:dimY
            k = index(i,j);
            A(k,k) = 1;
        end

        b(1:dimY) = TWest + 273.15;
        
    case 'Neumann'
        j = 1;
        
        for i = 1:dimY
            k = index(i,j);
            f = -lam(i,j)/(2*dx);
            A(k,k) = A(k,k) - 3*f;
            A(k,k+dimY) = A(k,k+dimY) + 4*f;
            A(k,k+2*dimY) = A(k,k+2*dimY) - 1*f;
        end

        b(1:dimY) = b(1:dimY) + beta;
        
    case 'Robin'
        j = 1;
        
        for i = 1:dimY
           k = index(i,j); 
           f = lam(i,j)/(2*dx);
           A(k,k) = A(k,k) + (alpha+3*f);
           A(k,k+dimY) = A(k,k+dimY) - 4*f;
           A(k,k+2*dimY) = A(k,k+2*dimY) + 1*f;
        end
       
        b(1:dimY) = b(1:dimY) + alpha*Tinf;
        
end

 

%------------------------------------------------------------------------
% Constructing Matrix A;

 

% central differencing sceme


for j = 1:dimX
    for i = 1:dimY
        k = index(i,j);
        if A(k,:) == zeros(dimA,1)
            A(k,k) = -(lam(i,j-1)+lam(i,j+1))/dx^2-(lam(i-1,j)+lam(i+1,j))/dy^2;
            A(k,k+1) = lam(i+1,j)/dy^2;
            A(k,k-1) = lam(i-1,j)/dy^2;
            A(k,k+dimY) = lam(i,j+1)/dx^2;
            A(k,k-dimY) = lam(i,j-1)/dx^2;
        end
    end
end

        
%------------------------------------------------------------------------      
% Solving the linear system (use '\' operator)
        
T = A\b;

 

Tmat = reshape(T,[dimX,dimY]);
        
%------------------------------------------------------------------------
% Ploting Results

 

% plot meshgrid (optional)
% figure(1)
 
% Do a surface plot
figure(1)
[X,Y] = meshgrid(0:dx:l,h:-dy:0);
figure(1);
surf(X,Y, Tmat);
title('2D Steady Heat Transfer');
colormap(jet);
colorbar;
xlabel('length l [m]');
ylabel('height h [m]');
zlabel('temperature T [K]');




% Do a contour plot
figure(2) 
contour(X,Y,Tmat);
title('2D Steady Heat Transfer');
colormap(jet);
colorbar;
xlabel('length l [m]');
ylabel('height h [m]');
zlabel('temperature T [K]')


 
%Do a color plot
figure(3)
pcolor(X,Y,Tmat);
title('2D Steady Heat Transfer');
colormap(jet);
colorbar;
xlabel('length l [m]');
ylabel('height h [m]');
zlabel('temperature T [K]');
