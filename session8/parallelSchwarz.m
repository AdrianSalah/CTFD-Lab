clear all;
close all;
clc;

%% Init Parameters

N = 13;         % number of points
dx = 1/(N-1);   % spatial resolution
delta = 1;      % normalized overlap (overlap in units of one dx)
Np = [2 4 6 12];         % number of subdomains
sol = zeros(N,1);   % solution vector

for i = 1:length(Np)
    % safe solution of each subdomain in a container
    u = cell(1, Np(i));
    A = cell(1, Np(i));
    b = cell(1, Np(i));

    %% Init problem for each subdomain and compute init solution (for each subdomain)

    for j=1:Np(i)

       if(j==1)
          A{j} = PoissonInit( (N-1)/Np(i)+delta+1, dx );
          b{j} = PoissonRHS( (N-1)/Np(i)+delta+1, dx, 0);
          u{j} = A{j}\b{j}';

       elseif (j==Np(i))
           A{j} = PoissonInit( (N-1)/Np(i)+delta+1, dx );
           b{j} = PoissonRHS( (N-1)/Np(i)+delta+1, dx, (j-1)/Np(i)-delta*dx );
           u{j} = A{j}\b{j}';

       else
           A{j} = PoissonInit( (N-1)/Np(i) + 2*delta+1, dx );
           b{j} = PoissonRHS( (N-1)/Np(i) + 2*delta+1, dx, (j-1)/Np(i)-delta*dx );
           u{j} = A{j}\b{j}';
       end

    end

    u0 = u;   


    %% Execute Parallel Schwarz Algorithm

    res = Inf;
    n = zeros(size(Np));
    maxtimestep = 1000;
    convergence_limit = 1.0e-6;

    Afull = PoissonInit(N,dx);
    bfull = PoissonRHS(N,dx,0);
    
    while (res > convergence_limit && n(i) < maxtimestep)


        for j=1:Np(i)

            b{j} = ExchangeBC(j, u0, b{j}, Np(i), delta);
            % spmd
            u{j} = A{j}\b{j}';
            % end
        end

        sol = combineSolution(u,N,Np(i));

        res = norm(bfull'-Afull*sol');

        u0 = u; 
        n(i) = n(i) + 1;
    end
    sprintf('Number of iterations for %d subdomains: %d',Np(i), n(i))
    
    %% plot solution
    x = 0:dx:1;     % spacial arry
    plot(x, sol);
    title('Temperature distribution 1D Poisson equation')
    xlabel('x')
    ylabel('T')

end