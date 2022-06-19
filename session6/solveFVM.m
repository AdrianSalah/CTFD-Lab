function T = solveFVM(T, X, Y, dimX, dimY, boundary, TD, alpha, Tinf, dt, tend, timeintegrationType, theta, simulationType)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File solveFVM.m
%
% This routine set up the linear system and solve it for every timestep
%
% input
% T        Spatial Matrix T
% X         Matrix x coordinates 
% Y         Matrix y coordinates
% boundary  String vector. Boundary types.
% TD        Temperature for each boundary (if Dirichlet)
% alpha     convective heat transfer coefficient
% Tinf      Temperature of the surrouding fluid 
% dt        timestepsize
% tend      end of time
% theta     parameter for theta scheme
%
% output
% T         Temperature field (matrix) for steady case
%           OR array of matrices 
%           Temperature field at every timestep) for transient case

% Index maps the node position to the correct linear equation
index = @(ii, jj) ii + (jj-1) * size(T, 1);


% B is the right-hand side of the linear system
B = zeros(size(T));

% create initial temperature field
if strcmp(boundary.east, 'Dirichlet')
     T(:, end) = TD.east;
end

if strcmp(boundary.west, 'Dirichlet')
     T(:, 1) = TD.west;
end

if strcmp(boundary.south, 'Dirichlet')
     T(end, :) = TD.south;
end

if strcmp(boundary.north, 'Dirichlet')
     T(1, :) = TD.north;
end


% set up the system matrix A
A = zeros(size(T, 1) * size(T, 2));

for i = 1:size(T, 1)
    for j = 1:size(T, 2)
        % Fill the system matrix and the right-hand side for node (i,j)
        [stecil,b] = stamp(i, j, X, Y, TD, alpha, Tinf, boundary);
        
        A(index(i,j), :) = stecil;
        B(i, j) = b;
    end
end

switch simulationType
    case 'steady'
%        spy(A);
        S = sparse(A);
%         
%         whos S
%         
%         whos A
%         tic
         T(:) = S\B(:);
%         toc
%         
%         tic
%         T(:) = A\B(:);
%         toc
%        b = B(:);
       % arrays for plotting res vs. iteration
%         resVec = [];
%         itVec = [];
       
       
        %JACOBI iteration
%         D = diag(A);
%        tolerance = 1.0e-6;
%        iterations = 100;
%        x = zeros(size(A,1),1);
%        % initial residual
%        res = b - A*x;
%        it = 0;
%        Dinv = 1./D;
%        % main Jacobi loop
%        while((norm(res)/norm(b) > tolerance) && (it < iterations))
%         x = x + diag(Dinv)*(b - A*x);
%         res = b - A*x;
% 
%         it = it + 1;
%         
%         % append current iteration and residual for plotting 
%         itVec = [itVec it];
%         resVec = [resVec norm(res)];
%         
%        end

%        tolerance = 1.0e-6;
%        iterations = 100;
%        relaxation = 1.9;
% 
% D = diag(A);
% E = tril(A,-1);
% F = triu(A,1);
% dimA = size(A,1);
% Dinv = 1./D;
% % initial guess
%        x = zeros(dimA,1);
%        % initial residual
%        res = b - A*x;
%        it = 0;
%        % main SOR loop
%        while(norm(res)/norm(b) > tolerance && it < iterations)
%            xnew = zeros(dimA,1);
%            for i = 1:dimA
%                xnew(i) = (1-relaxation)*x(i) + relaxation*Dinv(i)*( b(i) - E(i,:)*xnew - F(i,:)*x );
%            end
%            res = b - A*xnew;
%            x = xnew;  
%            it = it + 1;
%            
%        % append current iteration and residual for plotting 
%        itVec = [itVec it];
%        resVec = [resVec norm(res)];
%        end
% 
% 
%        figure(iterations)
%        plot(itVec,resVec, 'r-');
%        title('2D Heat Equation FVM, SOR Method \omega = 1.9')
%         xlabel('iteration')
%         ylabel('residual norm')
%         grid on
%         
%         hold off


    case 'transient'

        % get rhs vector  
        bvec = reshape(B, dimX*dimY, 1);

        % solve unsteady case
        switch timeintegrationType

%       NOT NEEDED!            
%             case 'explicitEuler'
%                 T = solveExplicitEuler(dt, tend, T, A, bvec, dimX, dimY);
% 
%             case 'implicitEuler'
%                 T = solveImplicitEuler(dt, tend, T, A, bvec, dimX, dimY);

            case 'theta'
                T = solveTheta(dt, tend, T, A, bvec, theta, dimX, dimY);

            case 'rungeKutta4'
                T = solve4RungeKutta(dt, tend, T, A, bvec, dimX, dimY);

            otherwise
                assert(false, 'false time integration specified: %s', timeintegrationType);
        end

end

