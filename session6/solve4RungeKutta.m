function [Tarr] = solve4RungeKutta(dt, tend, T0, A, b, dimX, dimY)
%SOLVEEXPLICIT Summary of this function goes here
%   returns:        
%   - Tarr: container of solution matrix at each timestep
%   parameters
%   - dt:             timestepsize
%   - tend:           end of time
%   - T0:             initial temperature distribution matrix
%   - A:              spatial discretization matrix  
%   - b:              rhs vector
%   - theta:          parameter for weighted average scheme
%   - dimX/ dimY:     dimension of solution matrix 


numTimesteps = tend/dt;
Tarr = zeros(dimX, dimY, numTimesteps);

% initial timestep
Tarr(:,:,1) = T0;



for i=1:numTimesteps-1
      Tvec_old = reshape(Tarr(:,:,i),dimX*dimY, 1);
      % calculate predictors and correctors
      T1 = Tvec_old + 0.5*dt*A*Tvec_old - b*dt/2;
      T2 = Tvec_old + 0.5*dt*T1 - b*dt/2;
      T3 = Tvec_old + dt*A*T2 - b*dt;
      Tvec_new = Tvec_old + 1/6*dt*(A*Tvec_old + 2*A*T1 + 2*A*T2 + A*T3) - b*dt;
      % reshape T vector of new timestep into matrix and store in array
      Tarr(:,:,i+1) = reshape(Tvec_new, dimX, dimY);   
end

