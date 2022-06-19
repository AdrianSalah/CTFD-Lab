function [T] = Lin_solver(A, b, Tinit, method, tolerance, relaxation,iterations)
%The function (that replaces the line T = A\B) you have to create reads:
% Inputs:
% where A is a NxN Matrix
% b is a Nx1 vector
% Tinit is a Nx1 vector (initial solution or guess)
% method is 'JAC', 'GAU', 'SOR'
% 
% tolerance (defined in Eq. 5.10) ... I will consider 1e-3 and 1e-4
% relaxation (for SOR)
% 
% iterations (number of maximum iterations)


switch method
    case 'JAC'
        D = diag(A);
        Dinv = 1./D;
        T = Tinit;
        res = b - A*T;
        it = 0;
        while(norm(res)/norm(b) > tolerance && it < iterations)
            T =  T + diag(Dinv)*(b - A*T);
            res = b - A*T;
            it = it + 1;
        end
        
    case 'GAU' 
        D = diag(A);
        Dinv = 1./D;
        E = tril(A,-1);
        F = triu(A,1);
        dimA = size(A,1);
        
        T = Tinit;
        res = b - A*T;
        it = 0;
       
       while(norm(res)/norm(b) > tolerance && it < iterations)
           
           Tnew = zeros(dimA,1);
           
           for i = 1:dimA
               Tnew(i) = Dinv(i)*( b(i) - E(i,:)*Tnew - F(i,:)*T );
           end
           
           res = b - A*Tnew;
           T = Tnew;  
           it = it + 1;
       end
        
    case 'SOR'
        D = diag(A);
        Dinv = 1./D;
        E = tril(A,-1);
        F = triu(A,1);
        dimA = size(A,1);
        
        T = Tinit;
        res = b - A*T;
        it = 0;
       
       while(norm(res)/norm(b) > tolerance && it < iterations)
           
           Tnew = zeros(dimA,1);
           
           for i = 1:dimA
               Tnew(i) = (1-relaxation)*T(i) + relaxation*Dinv(i)*( b(i) - E(i,:)*Tnew - F(i,:)*T );
           end
           
           res = b - A*Tnew;
           T = Tnew;  
           it = it + 1;
       end
        
end



end