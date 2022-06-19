function bool = isWaveStable(dimY, dimX, h, l, c, dt)
%WAVSTABILITY Summary of this function goes here
%   Detailed explanation goes here


    % TODO: extend for different formfunctions!
    dy = h/dimY;
    dx = l/dimX;
    
    
    if dt < 1/c* ( 1/dx^2 + 1/dy^2)^(-1/2)
       bool = true;
       
    else
       bool = false;
    end
    
    
end

