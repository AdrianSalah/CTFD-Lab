clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to solve the 2D steady heat equation in a non-Cartesian Grid by
% the Finite Volumes Method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize variables

InitFVM

% initialize spatial Matrix T

T = zeros(dimY, dimX);


% set up the mesh

[X, Y] = setUpMesh(T, l, formfunction);

% Fill matrix A and vector B. Solve the linear system.

[T, T_vec] = solveFVM(T, X, Y, boundary, TD, alpha, Tinf, betaeast, betanorth, betawest,k, tempdev, time, solver);

% Make some plots
switch tempdev
    case 'steady'
        VisTemperature(formfunction, X, Y, T_vec, l)
    case 'transient'
        VisTemperature(formfunction, X, Y, T_vec,l)
       
end
    








