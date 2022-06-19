clear all;
close all;

%%  CTFD 1. Introduction
% Form Adrian Salah & Daniel Eckert

%Analytical solution for 1D Convection/Diffusion/Convection & Diffusion
%equations


%% Initialisation

% Manipulatable values
D               = 0.1;              % Diffusivity
L               = 50;               % Period
d               = -20:0.1:20;       % Domain length in m
u               = 1;                % velocity in m/s
t_end           = 10;               % End of time domain in s
delta_t         = 0.1;              % timestepsize in s
N               = 150;              % number of terms in fourier series, 
                                    % for higher N we get better approximation of the exact analytical solution
num_timesteps   = t_end/delta_t;

% Adapted values
m       = 1:1:N;       


%% Coumputation of the Fourier coefficients 
% define symbolic variable for integration
syms x

% calculate sub-integrals for ranges of picewise defined function

% 0<x<=1
a_1 = x*sin(m*pi()*x/L);
a_m1 = double(int(a_1,0,1));

% 1<x<=2
a_2 = (2-x)*sin(m*pi()*x/L);
a_m2 = double(int(a_2,1,2));

% do not perform calculation explicitly to safe computational time
% 2<x<=L
% a_m3 = m*0;

% sum up sub-integrals to obtain final array of the N fourier coefficients
a_m = 2/L*(a_m1+a_m2);

%% Convection

% approximation of function for convective case as finite fourier series
theta_tilde_convection = @(x, t) sum(a_m.*sin(pi/L*m.*(x-u*t)));

% initialize result array
resultArrayConvection = zeros((num_timesteps + 1),length(d));

% compute result for given time- and x-range
for tn = 0:num_timesteps
    t = tn*delta_t;
    resultArrayConvection(tn+1, :) = arrayfun(@ (x) theta_tilde_convection(x, t), d);
end

% plot function and export as gif-file
plotAnalytical(resultArrayConvection, d, t_end, delta_t, 'Convection');



%% Diffusion

% approximation of function for diffusive case as finite fourier series
theta_tilde_diffusion = @(x, t) sum(a_m.*exp(-D*(m*pi/L).^2*t).*sin(pi/L*m.*x));

% initialize result array
resultArrayDiffusion = zeros((t_end/delta_t+1),length(d));

% compute result for given time- and x-range
for tn = 0:num_timesteps
    t = tn*delta_t;
    resultArrayDiffusion(tn+1, :) = arrayfun(@ (x) theta_tilde_diffusion(x, t), d);
end

% plot function and export as gif-file
plotAnalytical(resultArrayDiffusion, d, t_end, delta_t, 'Diffusion');



%% Convection and Diffusion

% approximation of function for diffusive and convetive case as finite fourier series
theta_tilde_convdiff = @(x,t) sum(a_m.*exp(-D*(m*pi/L).^2*t).*sin(pi/L*m.*(x-u*t)));

% initialize result array
resultArrayConvDiff = zeros((t_end/delta_t+1),length(d));

% compute result for given time- and x-range
for tn = 0:num_timesteps
    t = tn*delta_t;
    resultArrayConvDiff(tn+1, :) = arrayfun(@ (x) theta_tilde_convdiff(x, t), d);
end

% plot function and export as gif-file
plotAnalytical(resultArrayConvDiff, d, t_end, delta_t, 'Convection_Diffusion');

