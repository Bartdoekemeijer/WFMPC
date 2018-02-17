% This script simulates a MPC using the linear model and
% applies the input to the nonlinear model.
% The axial induction of the first turbine changes and the objective is to
% regulate the change in rotor velocity of the downwind tubines to zero using
% the axial induction factor of these turbines.

clear; clc; close all;

addpath bin
addpath(genpath('WFSim'))

%% Initialize script
Wp.name      = '6turb_adm_turb';    % Choose which scenario to simulate. See 'bin/core/meshing.m' for the full list.

% Model settings (recommended: leave default)
scriptOptions.Projection        = 1;        % Solve WFSim by projecting away the continuity equation (bool). Default: false.
scriptOptions.Linearversion     = 1;        % Calculate linear system matrices of WFSim (bool).              Default: false.
scriptOptions.exportLinearSol   = 1;        % Calculate linear solution of WFSim (bool).                     Default: false.
scriptOptions.Derivatives       = 0;        % Compute derivatives, useful for predictive control (bool).     Default: false.
scriptOptions.exportPressures   = ~scriptOptions.Projection;   % Calculate pressure fields. Default: '~scriptOptions.Projection'

% Convergence settings (recommended: leave default)
scriptOptions.conv_eps          = 1e-6;     % Convergence threshold. Default: 1e-6.
scriptOptions.max_it_dyn        = 1;        % Maximum number of iterations for k > 1. Default: 1.

% Display and visualization settings
scriptOptions.printProgress     = 1;    % Print progress in cmd window every timestep. Default: true.
scriptOptions.printConvergence  = 0;    % Print convergence values every timestep.     Default: false.
scriptOptions.Animate           = 100;   % Plot flow fields every [X] iterations (0: no plots). Default: 10.
scriptOptions.plotMesh          = 0;    % Plot mesh, turbine locations, and print grid offset values. Default: false.

scriptOptions.Control           = 0;

% WFSim general initialization script
[Wp,sol,sys] ...
    = InitWFSim(Wp,scriptOptions);

% Initialize variables and figure specific to this script
sol_array = {}; % Create empty array to save 'sol' to at each time instant
CPUTime   = zeros(Wp.sim.NN,1); % Create empty matrix to save CPU timings
if scriptOptions.Animate > 0 % Create empty figure if Animation is on
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
end
if Wp.sim.startUniform==1
    scriptOptions.max_it = 1;               % Maximum n.o. of iterations for k == 1, when startUniform = 1.
else
    scriptOptions.max_it = 50;              % Maximum n.o. of iterations for k == 1, when startUniform = 0.
end

% Controller variables
global uc;
perturbatie      = 1*-.4;                      % Perturbation of beta1 (disturbance rejection)
controller       = struct;

% Set control signals constant
for kk=1:Wp.sim.NN
    Wp.turbine.input(kk).CT_prime   = 2*ones(Wp.turbine.N,1);
    Wp.turbine.input(kk).dCT_prime  = zeros(Wp.turbine.N,1);
    Wp.turbine.input(kk).beta       = 0.25*Wp.turbine.input(kk).CT_prime;
    Wp.turbine.input(kk).dbeta      = zeros(Wp.turbine.N,1);
end

% Performing forward time propagations
disp(['Performing ' num2str(Wp.sim.NN) ' forward simulations..']);
while sol.k < Wp.sim.NN
    
    tic;                    % Start stopwatch
    [sol,sys]      = WFSim_timestepping(sol,sys,Wp,scriptOptions); % forward timestep: x_k+1 = f(x_k)
    
    % Save sol to cell array
    sol_array{sol.k} = sol;
    
    if sol.k>2
        
        controller = Computecontrolsignal(Wp,sys,controller,sol,perturbatie,scriptOptions);
        
        for kk=1:Wp.turbine.N-1
            Wp.turbine.input(sol.k+1).beta(kk+1)     = Wp.turbine.input(1).beta(kk+1) + uc(kk);
            Wp.turbine.input(sol.k+1).CT_prime(kk+1) = Wp.turbine.input(1).CT_prime(kk+1) + 4*uc(kk);
        end
        
        Wp.turbine.input(sol.k+1).beta(1)            = Wp.turbine.input(1).beta(1) + perturbatie;
        Wp.turbine.input(sol.k+1).CT_prime(1)        = Wp.turbine.input(1).CT_prime(1) + 4*perturbatie;
        
    end
    CPUTime(sol.k) = toc;   % Stop stopwatch
    
    
    % Print progress, if necessary
    if scriptOptions.printProgress
        disp(['Simulated t(' num2str(sol.k) ') = ' num2str(sol.time) ...
            ' s. CPU: ' num2str(CPUTime(sol.k)*1e3,3) ' ms.']);
    end;
    
    % Plot animations, if necessary
    if scriptOptions.Animate > 0
        if ~rem(sol.k,scriptOptions.Animate)
            hfig = WFSim_animation(Wp,sol,hfig);
        end;
    end;
end;
disp(['Completed ' num2str(Wp.sim.NN) ' forward simulations. Average CPU time: ' num2str(mean(CPUTime)*10^3,3) ' ms.']);

controller.znl = [zeros(1,Wp.sim.NN); controller.znl];
for kk=1:Wp.sim.NN
    CT_prime(:,kk) =  Wp.turbine.input(kk).CT_prime;
end

seq   = [2 4 6 1 3 5];

figure(1);clf
ll = 0;
for kk=seq
    ll = ll + 1;
    subplot(2,3,ll)
    if kk~=1
        stairs(Wp.sim.time(1:end-1),controller.znl(kk,:),'b');hold on; %This is what you really measure
        str = strcat('$U_',num2str(kk,'%.0f'),'$');
        ylabel([str,' [m/s]'],'interpreter','latex');
        xlabel('$t$ [s]','interpreter','latex')
        grid;xlim([0 Wp.sim.time(end)]);
    end
    if ll==2;title({'Averaged rotor flow velocity $U_{i}$';' '},'interpreter','latex');end;
    if ll==1;annotation(gcf,'arrow',[0.017 0.08],[0.51 0.51]);end;
end

figure(2);clf
ll = 0;
for kk=seq
    ll = ll + 1;
    subplot(2,3,ll)
    stairs(Wp.sim.time(1:end-1),CT_prime(kk,:),'b');hold on; %This is what you really measure
    str = strcat('$CT''_',num2str(kk,'%.0f'),'$');
    ylabel([str,' [-]'],'interpreter','latex');
    xlabel('$t$ [s]','interpreter','latex')
    grid;xlim([0 Wp.sim.time(end)]);
    if ll==2;title({'Control signal $CT''_{i}$';' '},'interpreter','latex');end;
    if ll==1;annotation(gcf,'arrow',[0.017 0.08],[0.51 0.51]);end;
end



