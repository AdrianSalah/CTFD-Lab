function T = solveFVM(T, X, Y, boundary, TD)


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
%
% output
% T         Temperature field

% Index maps the node position to the correct linear equation

index = @(ii, jj) ii + (jj-1) * size(T, 1);


% B is the right-hand side of the linear system
B = zeros(size(T));
if strcmp(boundary.north, 'Dirichlet')
    B(1, :) = TD.north;
end
if strcmp(boundary.south, 'Dirichlet')
    B(end, :) = TD.south;
end
if strcmp(boundary.east, 'Dirichlet')
    B(:, end) = TD.east;
end
if strcmp(boundary.west, 'Dirichlet')
    B(:, 1) = TD.west;
end


% set up the system matrix A
A = zeros(size(T, 1) * size(T, 2));

for i = 1:size(T, 1)
    for j = 1:size(T, 2)
        % Fill the system matrix and the right-hand side for node (i,j)
        [A(index(i,j), :)] = ...
            stamp(i, j, X, Y);
    end
end


% solve the linear system
T(:) = A\B(:);




