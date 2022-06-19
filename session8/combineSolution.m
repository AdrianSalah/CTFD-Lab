function [sol] = combineSolution(u, N, Np)
%COMBINESOLUTION combines solution of all workers (subdomains) in final solution array
%   Parameters:
%   u:  container of solutions of subdomains
%   N:  number of overall points in domain
%   Np: number of workers(subdomains)

%   Returns:
%   sols:   final solution array

% number of elements taken from the N-1 subarrays for the final
% solution
nsub = (N-1)/Np;

% accuumulate first subdomain
sol(1:nsub) = u{1}(1:end-2);

% accumulate inner subdomains
for j = 2:(Np-1)
    sol((j-1)*nsub+1:j*nsub) = u{j}(2:end-2);
end

% accuumulate last subdomain (here we include one additional point)
sol((Np-1)*nsub+1:Np*nsub+1) = u{Np}(2:end);

end

