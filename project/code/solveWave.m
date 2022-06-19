function [p_curr, p_arr] = solveWave(p_curr, A, B, dimX, dimY, X, Y, TD, alpha, Tinf, boundary, tend, dt, phase_vel, dummyNodes, maxIndexBoundaryWest, IndexNorthBoundary)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File solveWave.m
%
% This routine solves the 2d wave equation at every timestep
%
% input
% p_curr        Spatial Matrix p_curr (current pressure field)
% X             Matrix x coordinates 
% Y             Matrix y coordinates
% boundary      String vector. Boundary types.
% TD            Temperature for each boundary (if Dirichlet)
% alpha         convective heat transfer coefficient
% Tinf          Temperature of the surrouding fluid 
% dt            timestepsize
% tend          end of time
% phase_vel     phase velocity of the wave
% dummyNodes    nodes which don't lie in the specified geometry     
%
% output
% p_curr         pressure field matrix of current timestep
% p_arr          container of pressure field matrices for every timestep


index = @(ii, jj) ii + (jj-1) * size(p_curr, 1);
% container of matrices (pressure field at every timestep) for transient case
p_arr = zeros(dimY, dimX, tend/dt); 
% initialize solution container for every timestep
p_arr(:,:,1) = p_curr; 

for i = 1:size(p_curr, 1)
    for j = 1:size(p_curr, 2)
        % ignore dummy nodes
        if dummyNodes(index(i,j))==1
            continue;
        else   
            % Fill the system matrix and the right-hand side for node (i,j)
            A(index(i,j),:) = setupSystemMatrix(i, j, X, Y, TD, alpha, Tinf, boundary, IndexNorthBoundary);
        end    
    end
end

p_past = p_curr;
p_future = zeros(size(p_curr));
    
phase_vel = 1;

 for timestep=1:tend/dt   
    
    for i = 1:size(p_curr, 1)
        for j = 1:size(p_curr, 2)
            
            if dummyNodes(index(i,j))==1
                continue;
            else
                dTdt_p = (p_curr(i,j) - p_past(i,j))/dt;

                % Fill the system matrix and the right-hand side for node (i,j)
                B(i,j) = setupRHS(i, j, X, Y, TD, alpha, Tinf, boundary, timestep, phase_vel, dTdt_p, IndexNorthBoundary);
                
            end
            
        end
        
    end
   
    % update variables
    p_future(:) = phase_vel^2*dt^2*(A*p_curr(:) - B(:)) + 2*p_curr(:) - p_past(:); 
    % store current timestep in solution array
    p_arr (:,:,timestep) = p_curr;
   
    p_past = p_curr;
    p_curr = p_future;
    p_curr = updateBC(dimX, dimY, p_curr, TD, boundary, timestep, maxIndexBoundaryWest,IndexNorthBoundary);
    
    

 end
end
