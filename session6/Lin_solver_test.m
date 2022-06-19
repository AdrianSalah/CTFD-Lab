% dimension of A
n = 100;
tolerance = 1.0e-6;
iterations = 100;
relaxation = 1.9;
b = rand(n,1);

% arrays for plotting res vs. iteration
resVec = [];
itVec = [];

%% create random diagonal dominant matrix 
fprintf('random diagonal dominant matrix \n');

A = rand(n);
% A should be diagonal dominant
A(logical(eye(size(A)))) = n*n; % A diagonal dominant => GS and J converge, SOR convergence not guaranteed

%% create random three-diagonal dominant matrix
fprintf('random three-diagonal dominant matrix \n');

A = diag(rand(1,n-1),1) + diag(rand(1,n-1),-1) + diag(n*n*ones(1,n)); % A diagonal dominant => GS and J converge, SOR convergence not guaranteed

%% create random full matrix (convergence not guaranteed)
fprintf('random full matrix');
A = rand(n);

%% for testing

% compute solution for test
sol = A\b;


% check if A is singular

if det(A) ~= 0
    fprintf('A not sigular :) \n') 
else
    error('A is singular! :( \n')
end

%% Jacobi Method

D = diag(A);


% check if Diagonal matrix is invertable
if all(abs(D) > 0) 
   fprintf('D is diagonizable :) \n') 

   Dinv = 1./D;
   
   % check if spectral radius is smaller than 1 
   if (max(eig((diag(Dinv)*(A-diag(D))))) < 1)
       fprintf('spectral radius smaller than 1 :) \n')
       % initial guess
        x = zeros(size(A,1),1);
       % initial residual
       res = b - A*x;
       it = 0;
       % main Jacobi loop
       while((norm(res)/norm(b) > tolerance) && (it < iterations))
        x = x + diag(Dinv)*(b - A*x);
        res = b - A*x;

        it = it + 1;
        
        % append current iteration and residual for plotting 
        itVec = [itVec it];
        resVec = [resVec norm(res)];
        
       end
       
   else
       error('spectral radius is bigger than 1! :(')
   end
          
else
    error('D not diagonizable');
    
end

%% SOR Methode (Gauß-Seidel = SOR for relaxation = 1)

D = diag(A);
E = tril(A,-1);
F = triu(A,1);
dimA = size(A,1); 

if all(abs(D) > 0) 
   fprintf('D is diagonizable \n') 

   Dinv = 1./D;
   if (max(eig(inv(D+E))) < 1)
   % initial guess
       x = zeros(dimA,1);
       % initial residual
       res = b - A*x;
       it = 0;
       % main SOR loop
       while(norm(res)/norm(b) > tolerance && it < iterations)
           xnew = zeros(dimA,1);
           for i = 1:dimA
               xnew(i) = (1-relaxation)*x(i) + relaxation*Dinv(i)*( b(i) - E(i,:)*xnew - F(i,:)*x );
           end
           res = b - A*xnew;
           x = xnew;  
           it = it + 1;
           
       % append current iteration and residual for plotting 
       itVec = [itVec it];
       resVec = [resVec norm(res)];
       end
       
   else
      error('spectral radius is bigger than 1! :(') 
   end
   
   
   
else
    error('D not diagonizable');
end

%% Test result of iterative method


if norm(x-sol) < 1.0e-6
    fprintf('correct solution :) \n')
    
else
    error('wrong result for iterative method!')
    
end
   

%% plot residual vs. iteration

plot(itVec,resVec, 'r-');
%title('Random diag dominant Matrix, Jacobi Method')
%title('Random diag dominant Matrix, GS Method')
%title('Random diag dominant Matrix, SOR Method, \omega = 1.9')
%title('Random three-diag dominant Matrix, Jacobi Method')
%title('Random three-diag dominant Matrix, GS Method')
%title('Random three-diag dominant Matrix, SOR Method \omega = 1.9')
title('Random full Matrix, Jacobi Method')
xlabel('iteration')
ylabel('residual norm')
grid on




