% This script simulates a MPC using the linear model and
% applies the input to the nonlinear model.
% The axial induction of the first turbine changes and the objective is to
% regulate the change in rotor velocity of the downwind tubines to zero using
% the axial induction factor of these turbines.

clear; clc; close all;

addpath bin

%% Initialize script
options.Projection    = 1;                      % Use projection (true/false). For MPC put it true.
options.Linearversion = 1;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 1;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

Wp.name             = 'ThreeTurbine_mpc';       % Meshing name (see "\bin\core\meshing.m")
Wp.Turbulencemodel  = 'WFSim3';

global uc;

perturbatie   = 1*-.4;                         % Perturbation of beta1 (disturbance rejection)
Animate       = 0;                          % Show 2D flow fields every x iterations (0: no plots)
plotMesh      = 0;                          % Show meshing and turbine locations
conv_eps      = 1e-6;                       % Convergence threshold
max_it_dyn    = 1;                          % Maximum number of iterations for k > 1

if options.startUniform==1
    max_it = 1;
else
    max_it = 50;
end

% WFSim general initialization script
[Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] ...
    = InitWFSim(Wp,options,plotMesh);

controller    = struct;
if Wp.turbine.N==2
    [a1,a2]         = deal(zeros(1,Wp.sim.NN));
else
    [a1,a2,a3]      = deal(zeros(1,Wp.sim.NN));
end


if Animate > 0
    scrsz = get(0,'ScreenSize');
    hfig  = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
end

%% Loop
CPUTime = zeros(Wp.sim.NN-1,1);
for k=1:Wp.sim.NN-1
    tic
    
    it        = 0;
    eps       = 1e19;
    epss      = 1e20;
    
    while ( eps>conv_eps && it<max_it && eps<epss );
        it   = it+1;
        epss = eps;
        
        if k>1; max_it = max_it_dyn; end
        
        [sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k)] = ...
            Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options);                      % Create system matrices
        [sol,sys] = Computesol(sys,input{k},sol,k,it,options);                      % Compute solution
        [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);            % Map solution to field
        
        if k>2
            
            % These are only used for plotting
            a1(k)           = input{k}.beta(1)/(input{k}.beta(1)+1);
            a2(k)           = input{k}.beta(2)/(input{k}.beta(2)+1);
            if Wp.turbine.N==3
                a3(k)       = input{k}.beta(3)/(input{k}.beta(3)+1);
            end
            if k==3
                controller.CT    = CT(:,k);
                controller.Power = Power(:,k);
            end
            
            controller         = Computecontrolsignal(Wp,sys,controller,sol,input{k},k,perturbatie,options);
            
            input{k+1}.beta(1) = input{2}.beta(1) + perturbatie;
            input{k+1}.beta(2) = input{2}.beta(2) + uc(1);
            if Wp.turbine.N==3
                input{k+1}.beta(3) = input{2}.beta(3) + uc(2);
            end
        end
    end
    CPUTime(k) = toc;
    
    if Animate > 0
        if ~rem(k,Animate)
            Animation;
        end;
    end;
end;


%% Plot
if Wp.turbine.N==2
figure(2);clf
subplot(2,1,1)
plot(Wp.sim.time(3:end-2),controller.znl(3:end-1));
hold on;grid;hline(controller.ss,'k--');
title('$T_2$','interpreter','latex')
ylabel('$z^1$ [m/s]','interpreter','latex');
subplot(2,1,2)
plot(Wp.sim.time(3:end-2),a1(3:end-1));hold on;
plot(Wp.sim.time(3:end-2),a2(3:end-1),'r');grid;
ylim([min(min(a1),min(a2))-.1 max(max(a1),max(a2))+.05]);
ylabel('a');xlabel('Time [s]');
title('\color{blue} a^1, \color{red} a^2');
else
figure(2);clf
subplot(3,1,1)
plot(Wp.sim.time(3:end-2),controller.znl(1,3:end-1));
hold on;grid;xlim([0 350])
hline(controller.ss(1)+controller.r(1,1),'k--')
title('$T_2$','interpreter','latex')
ylabel('$z^1$ [m/s]','interpreter','latex');
subplot(3,1,2)
plot(Wp.sim.time(3:end-2),controller.znl(2,3:end-1));
hold on;grid;xlim([0 350])
hline(controller.ss(2)+controller.r(2,1),'k--')
title('$T_3$','interpreter','latex')
ylabel('$z^2$ [m/s]','interpreter','latex');
subplot(3,1,3)
plot(Wp.sim.time(3:end-2),a1(3:end-1));hold on;
plot(Wp.sim.time(3:end-2),a2(3:end-1),'k');
plot(Wp.sim.time(3:end-2),a3(3:end-1),'r');grid;xlim([0 300])
ylim([min(min(a1),min(a2))-.1 max(max(a1),max(a2))+.05]);
ylabel('$a$','interpreter','latex');xlabel('$k$ [s]','interpreter','latex');
title('$a^1$ (blue), $a^2$ (black), $a^3$ (red)','interpreter','latex'); 

figure(3);clf
subplot(3,1,1)
plot(Wp.sim.time(3:end-2),Power(1,3:end-1));
hold on;grid;xlim([0 350])
title('$T_1$','interpreter','latex')
ylabel('$P^1$ [W]','interpreter','latex');
subplot(3,1,2)
plot(Wp.sim.time(3:end-2),Power(2,3:end-1));
hold on;grid;xlim([0 350])
plot(Wp.sim.time(3:end-2),controller.Pr(1,3:end-1),'k--')
title('$T_2$','interpreter','latex')
ylabel('$P^2$ [W]','interpreter','latex');
subplot(3,1,3)
plot(Wp.sim.time(3:end-2),Power(3,3:end-1));
hold on;grid;xlim([0 350])
plot(Wp.sim.time(3:end-2),controller.Pr(2,3:end-1),'k--')
title('$T_3$','interpreter','latex')
ylabel('$P^3$ [W]','interpreter','latex');

end