clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to solve the 2D steady heat equation in a non-Cartesian Grid by
% the Finite Volumes Method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize variables

InitFVM

% set up the mesh
[X, Y, dummyNodes, IndexInletBoundary, IndexNorthBoundary] = setUpMesh(dimX, dimY, l, h, formfunction);

%% Solve wave equation

[p, p_arr] = solveFVM(X, Y, dimX, dimY, boundary, TD, alpha, Tinf, dt, tend, dummyNodes, c, IndexInletBoundary, IndexNorthBoundary);

% Store Workspace to file 
result_filename = ['results_', shape, '_dt', num2str(dt) ,'_tend', num2str(tend), '_dimX', num2str(dimX), '_dimY', num2str(dimY), '.mat'];
save(result_filename)


        
%% Time and spectral analysis

num_samples = length(p_arr);
sample_freq = 1/dt;
freq_range = (0:num_samples-1)*(sample_freq/num_samples);


% define arrays for probe values
probe_inlet_pipe = p_arr(i_probe_inlet, j_probe_inlet, :);
probe_outlet_chamber = p_arr(i_probe_outlet, j_probe_outlet, :);

probe_inlet_pipe = reshape(probe_inlet_pipe, [1, num_samples]);
probe_outlet_chamber = reshape(probe_outlet_chamber, [1, num_samples]);



% plot amplitude over time inlet probe
figure(1)
plot(time(:), probe_inlet_pipe(:))
xlabel('t [s]')
ylabel('p [Pa]')
title(['Time plot inlet probe, ', shape, ' dt=', num2str(dt), ', tend=', num2str(tend)])
grid on

% plot power DFT over freq inlet probe
figure(2)
freq_probe_inlet_pipe = fft(probe_inlet_pipe);
PowerDFT_inlet = abs(freq_probe_inlet_pipe).^2/num_samples;
plot(freq_range, PowerDFT_inlet);
%plot(freq_range, abs(freq_probe_inlet_pipe)); 
xlabel('f [Hz]')
xlim([0,1])
ylabel('Power')
title(['Spectrum inlet probe, ', shape , ' dt=', num2str(dt), ', tend=', num2str(tend)])
grid on

% plot amplitude over time outlet probe
figure(3)
plot(time(:), probe_outlet_chamber(:))
xlabel('t [s]')
ylabel('p [Pa]')
title(['Time plot outlet probe, ', shape , ' dt=', num2str(dt), ', tend=', num2str(tend)])
grid on

% plot power DFT over freq outlet probe
figure(4)
freq_probe_outlet_chamber = fft(probe_outlet_chamber);
PowerDFT_outlet = abs(freq_probe_outlet_chamber).^2/num_samples;
%plot(freq_range,abs(freq_probe_outlet_chamber)); 
plot(freq_range, PowerDFT_outlet);
xlabel('f [Hz]')
xlim([0,1])
ylabel('Power')
title(['Spectrum outlet probe, ', shape, ' dt=', num2str(dt), ', tend=', num2str(tend)])

grid on


% determine time when signal leaves inlet tube
index_inlet_opening = floor(openingPos2*dimX);

for i=1:size(p_arr,3)   
   east_boundary = p_arr(:,index_inlet_opening,i);
   
   if any(abs(east_boundary(:)) > 1.0e-6)
      time_inlet_signal = i*dt;
      break
   end
end



% determine time when signal reaches outlet
for i=1:size(p_arr,3)   
   east_boundary = p_arr(:,dimX,i);
   
   if any(abs(east_boundary(:)) > 1.0e-6)
      time_outlet_signal = i*dt;
      break
   end
end

% determine mean, variance for inlet probe
mean_inlet_probe = mean(abs(probe_inlet_pipe));
max_inlet_probe = max(abs(probe_inlet_pipe));

% determine mean, variance for outlet probe 

mean_outlet_probe = mean(abs(probe_outlet_chamber));
max_outlet_probe = max(abs(probe_outlet_chamber));


%% plot mesh, geometry and probe positions

figure(5)
Z0 = 10*ones(dimY, dimX);
Z1 = 10*ones(dimY, dimX);
Z1(i_probe_inlet, j_probe_inlet) = -40;
Z1(i_probe_outlet, j_probe_outlet) = -20;
hold on
colormap('gray')
pcolor(X,Y,Z1)
pcolor(X,-Y,Z0)

plot([0 linspace(0, l, 50) l], [0 formfunction(linspace(0, 1, 50)) 0], 'k-','LineWidth',2)
plot([0 linspace(0, l, 50) l], [0 -formfunction(linspace(0, 1, 50)) 0], 'k-','LineWidth',2)
axis off;


%% Transient plots and videos 
% create transient videos
VisTemperature_transient(formfunction, X, Y, p_arr, l, h, dt, 'transient test')

%% Store Workspace to file 

filename = ['results_', shape, '_dt', num2str(dt) ,'_tend', num2str(tend), '_dimX', num2str(dimX), '_dimY', num2str(dimY), '.mat'];
save(filename)


%% Load Result file for Evaluation

%%%%%%%%%%%%%%% This is an exaple! 
%%%%%%%%%%%%%%% Please specify an already existing result and setup filename here

% clear
% setup_filename = 'setup_angledBackstep_dt0.1_tend10_dimX30_dimY20.mat'
% result_filename = 'result_angledBackstep_dt0.1_tend10_dimX30_dimY20.mat'
% load(result_filename)
% load(result_filename)
 





