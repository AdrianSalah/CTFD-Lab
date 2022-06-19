
function [b] = ExchangeBC(j, u0, b0, Np, delta)
%EXCHANGEBC returns
%   Detailed explanation goes here

b = b0;
if(j==1)
    b(end) = u0{j+1}(1+2*delta); % set the right BC for omega_1
    % b(1) = 0;
elseif(j==Np)
    b(1) = u0{j-1}(end-2*delta); % set left BC for omega_Np
    % b(end) = 0;
else
    b(1) = u0{j-1}(end-2*delta); % set left BC for omega_j
    b(end) = u0{j+1}(1+2*delta); % set right bC for omega_j
end

