function [p_curr, p_arr] = solveFVM(X, Y, dimX, dimY, boundary, TD, alpha, Tinf, dt, tend, phase_vel, dummyNodes, IndexInletBoundary, IndexNorthBoundary)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File solveFVM.m
%
% This routine set up the linear system and solve it for every timestep
%
% input
% p_curr    Spatial Matrix p_curr (current pressure field)
% X         Matrix x coordinates 
% Y         Matrix y coordinates
% boundary  String vector. Boundary types.
% TD        Temperature for each boundary (if Dirichlet)
% alpha     convective heat transfer coefficient
% Tinf      Temperature of the surrouding fluid 
% dt        timestepsize
% tend      end of time
% phase_vel phase velocity of wave
%
% output
% p_curr         pressure field matrix of current timestep
% p_arr          container of pressure field matrices for every timestep

% Index maps the node position to the correct linear equation
index = @(ii, jj) ii + (jj-1) * size(p_curr, 1);

p_curr = zeros(dimY, dimX);

% B is the right-hand side of the linear system
B = zeros(size(p_curr));

% create initial pressure field (with Dirichlet BC)
if strcmp(boundary.east, 'Dirichlet')
     p_curr(:, end) = TD.east;
end

if strcmp(boundary.west, 'Dirichlet')
    % west boundary is restricted to non-dummy nodes
    p_curr(IndexInletBoundary+1:dimY, 1) = TD.west(1);
end

if strcmp(boundary.south, 'Dirichlet')
     p_curr(end, :) = TD.south;
end

if strcmp(boundary.north, 'Dirichlet')
     p_curr(IndexNorthBoundary) = TD.north;
end


% set up the system matrix A
A = zeros(dimX * dimY);


% solve unsteady case with wave equation solver
[p_curr, p_arr] = solveWave(p_curr, A, B, dimX, dimY, X, Y, TD, alpha, Tinf, boundary, tend, dt, dummyNodes, phase_vel, IndexInletBoundary, IndexNorthBoundary);


end

