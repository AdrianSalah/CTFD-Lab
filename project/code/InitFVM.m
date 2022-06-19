
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the parameters of the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define type of shape function for domain (explenation below!)
shape = 'smoothOutlet'; % choose: 'backstep' or 'angledBackstep' or 'smoothOutlet' 

% l: length, h: height of domain
l = 10;
h = 5;

% Number of degrees of freedom (number of nodes per length)
dimX = 20;
dimY = 20;

% Parameter for Conjugated Heat Transfer (For Session 04)
alpha = 0;
Tinf = 90;

% Timestepsize and Endtime for unsteady case
dt = 0.1;
tend = 10;

% phase velocity of wave
c = 1;


time = dt:dt:tend; 

%%%%%% Not USED FOR PROJECT
% determines if excitation is over full simulation time or just till
% tend_excitation

% is_fulltime_excitation = false; 
% tend_excitation = 10;
% 
% 
% if is_fulltime_excitation
%     % sinus excitation for the full simulation time
%     time = dt:dt:tend;   
% else
%     % sinus excitaion for one period and dirichlet BC = 0 afterwards
%     time_excitation = dt:dt:tend_excitation;
%     time = [time_excitation  zeros(1,(tend-tend_excitation)/dt)];
% end


% Signal for west Dirichlet BC
in_signal = 2*rand(tend/dt, 1)-1; % random (white noise signal)
% in_signal = sin(2*pi*time);     % uniform sinus input


if ~isWaveStable(dimY, dimX, h, l, c, dt)
   error('CFL condition not fullfilled'); 
   %warning('CFL condition not fullfilled')
end

% specify Boundary conditions
boundary.south = 'Neumann';
boundary.north = 'Neumann';
boundary.east = 'Open';
boundary.west = 'Dirichlet';

% set boundary condtions
TD.north = 0;
TD.south = 0;
TD.west =  100*in_signal; 
TD.east = 0;


switch shape
    
%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT USED FOR PROJECT    
%     case 'straigthLine'
%         
%         formfunction = @(xnorm) h1;
%     
%     case 'linear'
% 
%         formfunction = @(xnorm) (1-xnorm)*h1/2 + xnorm*h2/2;
% 
%     case 'quadratic'
% 
%         c1 = h2+2*h1/2-2*hm;
%         c2 = 2*hm - 3*h1/2 - h2/2;
%         c3 = h1/2;
% 
%         formfunction = @(xnorm) c1*xnorm.^2 +c2*xnorm + c3;
        


    case 'backstep'
        % Inlet on the west side
        %           ___________
        %          | 
        %__________|
        %
        % Inlet          Outlet
        %__________
        %          |
        %          |___________

        bigRadius = 5;
        smallRadius = 2;
        openingPos2 = 0.3;
        syms y(x)
        y(x) = piecewise(x<openingPos2, smallRadius, x>=openingPos2, bigRadius);
        formfunction = @(xnorm) y(xnorm);

    case 'angledBackstep'
        % Inlet on the west side
        %             __________
        %            |
        %___________/ <- angledCorner
        %
        % Inlet          Outlet
        %___________
        %           \
        %            |__________

        bigRadius = 5;
        intermediateRadius = 3;
        smallRadius = 2;
        openingPos1 = 0.25;
        openingPos2 = 0.35;
        m = (intermediateRadius - smallRadius)/(openingPos2 - openingPos1);
        t = smallRadius - m*openingPos1;
        syms y(x)
        y(x) = piecewise(x<openingPos1, smallRadius, x>=openingPos2, bigRadius, m*x + t);
        formfunction = @(xnorm) y(xnorm);

    case 'smoothOutlet'
        % Inlet on the west side
        %           _____
        %          |
        %        __|
        % ____---  <-smoothCorner
        %
        % Inlet       Outlet
        % ____
        %     ---__
        %          |
        %          |_____

        bigRadius = 5;
        intermediateRadius = 3;
        smallRadius = 2;
        openingPos1 = 0.25;
        openingPos2 = 0.35;
        %c1 = intermediateRadius - smallRadius;
        %c2 = smallRadius ;  
        exponent = 5;
        a = (smallRadius - bigRadius)/(openingPos1^exponent - openingPos2^exponent);
        b = (smallRadius - a*openingPos1^exponent);
        syms y(x)
        %y(x) = piecewise(x<openingPos1, smallRadius, (x>=openingPos1) && (x < openingPos2), c1 * xnorm.^5 + c2, x >= openingPos2, bigRadius);
        y(x) = piecewise(x<openingPos1, smallRadius, x >= openingPos2, bigRadius,  a * x^exponent + b);
        formfunction = @(xnorm) y(xnorm);

end


% probes for spectral analysis 
j_probe_outlet = floor(0.75*dimX);
i_probe_outlet = floor(0.5*dimY);

j_probe_inlet = floor(0.2*dimX);
i_probe_inlet = floor((1-0.1)*dimY);

%% Store Setup to file 

setup_filename = ['setup_', shape, '_dt', num2str(dt) ,'_tend', num2str(tend), '_dimX', num2str(dimX), '_dimY', num2str(dimY), '.mat'];
save(setup_filename)



