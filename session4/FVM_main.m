clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to solve the 2D steady heat equation in a non-Cartesian Grid by
% the Finite Volumes Method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize variables

InitFVM

% initialize spatial Matrix T

T = zeros(dimY, dimX);

% set up the mesh

[X, Y] = setUpMesh(T, l, formfunction);

%% Fill matrix A and vector B. Solve the linear system.


T = solveFVM(T, X, Y, boundary, TD, alpha, Tinf);

%% Make some plots
VisTemperature(formfunction, X, Y, T, l)








