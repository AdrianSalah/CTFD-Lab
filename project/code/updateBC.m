function p_curr = updateBC(dimX, dimY, p_curr, TD, boundary,timestep, IndexInletBoundary, IndexNorthBoundary)
%UPDATEBC Summary of this function goes here
%   Detailed explanation goes here
    
    
    
    
    if strcmp(boundary.east, 'Dirichlet')
         p_curr(:, end) = TD.east;
    end

    if strcmp(boundary.west, 'Dirichlet')
        %    
        p_curr(IndexInletBoundary+1:dimY, 1) = TD.west(timestep);
    end

    if strcmp(boundary.south, 'Dirichlet')
         p_curr(end, :) = TD.south;
    end

    if strcmp(boundary.north, 'Dirichlet')
         p_curr(IndexNorthBoundary) = TD.north;
    end

end

