clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to solve the 2D steady heat equation in a non-Cartesian Grid by
% the Finite Volumes Method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize variables

InitFVM

% initialize spatial Matrix T
% initialize Temperature Matrix with Dirichlet BC!
T = zeros(dimY, dimX);





% set up the mesh

[X, Y] = setUpMesh(T, l, formfunction);

%% Fill matrix A and vector B. Solve the linear system.

% --- SOLVE ---

switch simulationType
    
    case 'steady'
        T = solveFVM(T,X,Y,dimX,dimY,boundary,TD,alpha,Tinf,dt,tend,timeintegrationType, theta, simulationType);
        
    case 'transient'
        Tarr = solveFVM(T,X,Y,dimX,dimY,boundary,TD,alpha,Tinf,dt,tend,timeintegrationType, theta, simulationType);
        
    otherwise
        error(false, 'false simulation type specified: %s', simulationType);
end






%% Make some plots
switch simulationType
    
    case 'steady'
        VisTemperature(formfunction, X, Y, T, l)
        
    case 'transient'
        % create videos for transient case
        VisTemperature_transient(formfunction, X, Y, Tarr, l, 'transient test')
end


%% Plot Storage and comp. time (sparse vs. full matrix)

nMatr           = [10, 100, 1000 10000];

tSparse         = [0.001229,0.001519, 0.007084, 0.074505];
storageSparse   = 0.001*[688, 12616, 147080, 1492936]; % MB

tFull           = [0.000220,0.000665, 0.060072, 13.042308];
storageFull     = 0.001*[648, 80000, 8388608, 800000000]; % MB



figure(1)
%plot storage

semilogx(nMatr,storageSparse, 'b--');
xlabel('n')
ylabel('storage [MB]')
ylim([1.0e-4 1.0e4]);
grid on
hold on
semilogx(nMatr,storageFull, 'r-');

hold off

figure(2)
% plot time
xlabel('n')
ylabel('t [sec]')
%plot(nMatr,tSparse)
semilogx(nMatr,tSparse, 'b--');
xlabel('n')
ylabel('t [sec]')
grid on
hold on
semilogx(nMatr,tFull, 'r-');

