function [Tarr] = solveImplicitEuler(dt, tend, T0, A, b, dimX, dimY)
%SOLVEEXPLICIT Summary of this function goes here
%   returns:        
%   - Tarr: container of solution matrix at each timestep
%   parameters
%   - dt:             timestepsize
%   - tend:           end of time
%   - T0:             initial temperature distribution matrix
%   - A:              spatial discretization matrix  
%   - b:              rhs vector
%   - dimX/ dimY:     dimension of solution matrix 


numTimesteps = tend/dt;
Tarr = zeros(dimX, dimY, numTimesteps);
% initial timestep
Tarr(:,:,1) = T0;

Astar = eye(dimX*dimY) - dt * A;

for i=1:numTimesteps-1
    Tvec_old = reshape(Tarr(:,:,i),dimX*dimY, 1);
    bstar = Tvec_old - b*dt;
    % solve linear system in each timestep
    Tvec_new = Astar\bstar;
    % reshape T vector of new timestep into matrix and store in array
    Tarr(:,:,i+1) = reshape(Tvec_new, dimX, dimY);
    
end

