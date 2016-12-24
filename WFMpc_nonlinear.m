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

%Wp.name       = 'TwoTurbine_mpc';       % Meshing name (see "\bin\core\meshing.m")
Wp.name       = 'ThreeTurbine_mpc';       % Meshing name (see "\bin\core\meshing.m")

global uc;

perturbatie   = -.4;                         % Perturbation of beta1
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
            %
            
            controller         = Computecontrolsignal(Wp,sys,controller,sol,input{k},k,perturbatie);
            
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
title('$\overline{U^r}$ of $T_2$','interpreter','latex')
ylabel('$\overline{U^r}$ [m/s]','interpreter','latex');
subplot(2,1,2)
plot(Wp.sim.time(3:end-2),a1(3:end-1));hold on;
plot(Wp.sim.time(3:end-2),a2(3:end-1),'r');grid;
ylim([min(min(a1),min(a2))-.1 max(max(a1),max(a2))+.05]);
ylabel('a');xlabel('Time [s]');
title('\color{blue} a_1, \color{red} a_2');
else
figure(2);clf
subplot(3,1,1)
plot(Wp.sim.time(3:end-2),controller.znl(1,3:end-1));
hold on;grid;xlim([0 700])
hline(controller.ss(1),'k--')
title('$\overline{U^r}$ of $T_2$','interpreter','latex')
ylabel('$\overline{U^r}$ [m/s]','interpreter','latex');
subplot(3,1,2)
plot(Wp.sim.time(3:end-2),controller.znl(2,3:end-1));
hold on;grid;xlim([0 700])
hline(controller.ss(2),'k--')
title('$\overline{U^r}$ of $T_3$','interpreter','latex')
ylabel('$\overline{U^r}$ [m/s]','interpreter','latex');
subplot(3,1,3)
plot(Wp.sim.time(3:end-2),a1(3:end-1));hold on;
plot(Wp.sim.time(3:end-2),a2(3:end-1),'k');
plot(Wp.sim.time(3:end-2),a3(3:end-1),'r');grid;xlim([0 700])
ylim([min(min(a1),min(a2))-.1 max(max(a1),max(a2))+.05]);
ylabel('$a$','interpreter','latex');xlabel('$t$ [s]','interpreter','latex');
title('$a_1$ (blue), $a_2$ (black), $a_3$ (red)','interpreter','latex');  
end