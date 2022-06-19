function [T , T_vec] = solveFVM(T, X, Y, boundary, TD, alpha, Tinf, betaeast, betanorth, betawest,k,tempdev, time, solver)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File solveFVM.m
%
% This routine set up the linear system and solve it
%
% input
% T         Spatial Matrix T
% X         Matrix x coordinates 
% Y         Matrix y coordinates
% boundary  String vector. Boundary types.
% TD        Temperature for each boundary (if Dirichlet)
% alpha     convective heat transfer coefficient
% Tinf      Temperature of the surrouding fluid 
% tempdev   contains information: steady/transient
% output
% T         Temperature field
% T_vec     Temperature for each time step
 
% Index maps the node position to the correct linear equation

index = @(ii, jj) ii + (jj-1) *  size(T, 1);


% B is the right-hand side of the linear system
B = zeros(size(T, 1) * size(T, 2),1);       

% set up the system matrix A
A = zeros(size(T, 1) * size(T, 2));
T_vec = zeros (size(T,1)*size(T,2),time(2)+1); %%allocate space for following timesteps

switch tempdev
    case 'steady'

for i = 1:size(T, 1)
    for j = 1:size(T, 2)
%        disp(index(i,j));
%        disp(index(j,i));
        % Fill the system matrix and the right-hand side for node (i,j)
        [temp1, temp2] = stamp(i, j, X, Y, TD, alpha, Tinf, boundary, betaeast, betanorth, betawest,k);
            [A(index(i,j), :)] = ... % fills line of Matrix A
                temp1;
            B(index(i,j)) = temp2;
        
    end
end
A = sparse(A);
T(:) = A\B(:);
    
    case 'transient'
        
if ~strcmp(boundary.west, 'Dirichlet') && ~strcmp(boundary.north, 'Dirichlet') && ~strcmp(boundary.east, 'Dirichlet') && ~strcmp(boundary.south, 'Dirichlet')

T (:,1) = 100; %% define default temperature distribution for no Dirichlet condition
else
% TW = zeros(size(T));
% TS = zeros(size(T));
% TN = zeros(size(T));
% TE = zeros(size(T));

if strcmp(boundary.west, 'Dirichlet')
    
    T(:,1) = TD.west;
end
if strcmp(boundary.north, 'Dirichlet')
    
    T(1,:) = TD.north;
end
if strcmp(boundary.east, 'Dirichlet')
    
    T(:,size(T,2)) = TD.east;
end
if strcmp(boundary.south, 'Dirichlet')
    
    T(size(T,1),:) = 20;
end

% T = (TW + TS +TN + TE);
% T(1,1) = 0.5 *T(1,1);
% T(size(T,1),1) = 0.5*T(size(T,1),1);
% T(1,size(T,2)) = 0.5*T(1, size(T,2));
% T(size(T,1), size(T,2)) = 0.5*T(size(T,1), size(T,2));
end


T_vec(:,1) = T(:);
for t = 1:time(2) 
    for i = 1:size(T, 1)
        for j = 1:size(T, 2)
        % Fill the system matrix and the right-hand side for node (i,j)
        [temp1, temp2] = stamp(i, j, X, Y, TD, alpha, Tinf, boundary, betaeast, betanorth, betawest,k);
            [A(index(i,j), :)] = ... % fills line of Matrix A
                temp1;
            B(index(i,j)) = temp2;
            
       end
    end
    
   switch solver.type
       case 'theta'
           theta = solver.type2;
           T(:) = sparse((eye(size(T,1)*size(T,2))*(1/time(3))-theta*A))\sparse((-B(:)+(1-theta)*A*T(:)+(1/time(3) *T(:))));
           T_vec (:,t+1) = T(:);
      
       case 'RK4'
           Td= T(:) + 0.5* time(3 )*A*T(:) -0.5*B(:)*time(3);
           Tdd = T(:) + 0.5*time(3)*Td -0.5*B(:)*time(3);
           Tddd= T(:) + time(3)*A*Tdd- B(:)*time(3);
           T(:) = T(:) + (1/6)*time(3)*(A*T(:)+2*A*Td+2*A*Tdd+A*Tddd)-B(:)*time(3);
           T_vec (:,t+1) = T(:);
   end
%          
end
end
% solve the linear system






