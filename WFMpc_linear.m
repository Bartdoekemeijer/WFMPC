%% Simulate the linear model with MPC

% Compute steady state solution with WFSim (Projection=1&&Linearversion=1)
% and save the workspace using: clear hfig;save('Data_WFSim\nonlinearsol').
% Then run this script to compute a MPC using the linear model and 
% apply this to the same model. 
% The axial induction of the first turbine changes and the objective is to
% regulate the change in rotor velocity of the second turbine to zero using 
% the axial induction factor of the second turbine.

clear; clc; close all;

Animate   = 0;
Movie     = 0;
Truncated = 0;  % Sparsified system matrix
Control   = 1;

nu      = Wp.Nu;
nv      = Wp.Nv;
np      = Wp.Np;
n1      = (Wp.mesh.Nx-3)+(Wp.mesh.Nx-2);    % # cells in x-direction for one y position
n2      = size(Wp.mesh.yline{1},2);         % # cells in y-direction where the turbine is
h       = Wp.sim.h;
L       = 250;
time    = (0:h:L);
NN      = length(time);

x       = zeros(size(sys.Qsp,2),NN);            % alpha state
x1      = zeros(size(sys.Qsp,2),NN);            % alpha state truncated sys
xnc     = zeros(size(sys.Qsp,2),NN);            % alpha state noncontrolled sys
xnc1    = zeros(size(sys.Qsp,2),NN);            % alpha state noncontrolled truncated sys
y       = zeros(nu+nv,NN);                      % [u;v]
y1      = zeros(nu+nv,NN);                      % [u;v] truncated sys
ync     = zeros(nu+nv,NN);                      % [u;v] uncontrolled sys
ync1    = zeros(nu+nv,NN);                      % [u;v] uncontrolled truncated sys
z       = zeros(1,NN);                          % Performance channel
z1      = zeros(1,NN);                          % Performance channel truncated sys
znc     = zeros(1,NN);                          % Performance channel uncontrolled sys
znc1    = zeros(1,NN);                          % Performance channel uncontrolled truncated sys
nz      = size(z,1);

% Control inputs
deltabeta   = 0.2;  % This is a disturbance on first turbine
deltayaw    = 0;
w           = [deltabeta*ones(1,NN);-0*deltabeta*ones(1,NN) ; ...
    deltayaw*ones(1,NN);-deltayaw*ones(1,NN)]; % Perturbation signals (note that w(2)=beta2 is the control variable)
r           = 0*ones(1,NN); % Reference to sum of rotor velocity of turbine 2
uc          = zeros(1,NN); % Control signal: beta2

% Flow fields
[du,dv,ul,vl]     = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny,NN));
[du1,dv1,ul1,vl1] = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny,NN));

% System matrices
A          = sys.Etl\sys.Atl;
Af         = A;
B          = sys.Etl\sys.Btl;
Bf         = B;

A(A<10e-2) = 0;             % Make system matrix sparse

Ctl        = sys.Qsp;
Cz         = zeros(1,nu+nv);
Cz((Wp.mesh.xline(2)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{2}(1)-1:...
    (Wp.mesh.xline(2)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{2}(end)-1) = 1;
C          = Cz*Ctl;
D          = zeros(size(C,1),size(w,1));

% ss(Atl,Btl,Ctl,0,Etl,h);      % Model from w to y
% ss(A,B,Ctl,0,h);              % Model from w to y1
% ss(A,B,C,0,h);                % Model from w to z

if Control>0
    if Truncated>0
        sysMpc    = ss(full(A),full(B),C,D,Wp.sim.h);  % Truncated
    else
        sysMpc    = ss(full(Af),full(Bf),C,D,Wp.sim.h);  % Full
    end
end

display(' ')
display('Compute mpcstate. Will take some time..')
display(' ')

if Control>0
    model   = setmpcsignals(sysMpc,'MD',[1 3 4],'MV',2,'MO',1);
    MV      = struct('Min',-.2,'Max',.2,'RateMin',-.01,'RateMax',.01);
    if Truncated>0
        mpcobj  = mpc(model,Wp.sim.h,10,1,[],MV);
    else
        mpcobj  = mpc(model,Wp.sim.h,25,1,[],MV);
    end
    xmpc    = mpcstate(mpcobj);
    options = mpcmoveopt;
end

% [L,M,A1,Cm1] = getEstimator(mpcobj); % Retrieve Kalman gain

% Figure
if Animate==1
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
    if Movie==1;
        ha=gcf;
        aviobj = VideoWriter(['..\..\Data_WFSim\Compare_linearsys_with_truncated','.avi']);
        aviobj.FrameRate=5;
        open(aviobj)
    end
end
display(' ')
display('Start the time-loop')
display(' ')

% Time loop
for kk=1:NN-1
    tic
    yaw_angles = .5*Wp.turbine.Drotor*exp(1i*w(Wp.turbine.N+1:end,kk)*pi/180);  % Yaw angles
    
    y(:,kk)     = Ctl*x(:,kk);      % Output full system
    y1(:,kk)    = Ctl*x1(:,kk);     % Output truncated system
    ync(:,kk)   = Ctl*xnc(:,kk);    % Output uncontrolled system
    ync1(:,kk)  = Ctl*xnc1(:,kk);   % Output uncontrolled truncated system
    
    z(:,kk)    = Cz*y(:,kk);    % Sum rotor velocities turbine 2
    z1(:,kk)   = Cz*y1(:,kk);   % Sum rotor velocities turbine 2 truncated
    znc(:,kk)  = Cz*ync(:,kk);  % Sum rotor velocities uncontrolled turbine 2
    znc1(:,kk) = Cz*ync1(:,kk); % Sum rotor velocities uncontrolled truncated turbine 2
    
    % Compute control action second turbine
    if Control>0
        if Truncated>0
            uc(kk) = mpcmove(mpcobj,xmpc,z1(:,kk),r(:,kk),[w(1,kk);w(3:4,kk)],options);
        else
            uc(kk) = mpcmove(mpcobj,xmpc,z(:,kk),r(:,kk),[w(1,kk);w(3:4,kk)],options);
        end
    end
    % Update plant states
    
    x(:,kk+1)   = sys.Etl\(sys.Atl*x(:,kk) + sys.Btl*[w(1,kk);uc(kk);w(3:4,kk)]);
    x1(:,kk+1)  = A*x1(:,kk) + B*[w(1,kk);uc(kk);w(3:4,kk)];
    xnc(:,kk+1) = sys.Etl\(sys.Atl*xnc(:,kk) + sys.Btl*w(:,kk));
    xnc1(:,kk+1)= A*xnc1(:,kk) + B*w(:,kk);
    
    du(3:end-1,2:end-1,kk)  = reshape(y(1:nu,kk),Wp.mesh.Ny-2,Wp.mesh.Nx-3)';
    dv(2:end-1,3:end-1,kk)  = reshape(y(nu+1:nu+nv,kk),Wp.mesh.Ny-3,Wp.mesh.Nx-2)';
    du1(3:end-1,2:end-1,kk) = reshape(y1(1:nu,kk),Wp.mesh.Ny-2,Wp.mesh.Nx-3)';
    dv1(2:end-1,3:end-1,kk) = reshape(y1(nu+1:nu+nv,kk),Wp.mesh.Ny-3,Wp.mesh.Nx-2)';
    
    ul(:,:,kk)  = sol.u + du(:,:,kk);
    vl(:,:,kk)  = sol.v + dv(:,:,kk);
    ul1(:,:,kk) = sol.u + du1(:,:,kk);
    vl1(:,:,kk) = sol.v + dv1(:,:,kk);
    
    %z(:,kk) - sum(du1(Wp.mesh.xline(2),Wp.mesh.yline{2},kk)) % Has to be zero
    
    toc
    
    if Animate==1
        if ~rem(kk,2)
            
            subplot(2,2,[1 3]);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',ul(:,:,kk),'Linecolor','none');  colormap(hot);
            caxis([min(min(ul(:,:,kk))) max(max(ul(:,:,kk)))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            text(0,Wp.mesh.ldxx2(end,end)+80,['Time ', num2str(time(kk),'%.1f'), 's']);
            xlabel('y [m]')
            ylabel('x [m]');
            title('u_l [m/s]');
            hold off;
            
            subplot(2,2,[2 4]);
            
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',ul1(:,:,kk),'Linecolor','none');  colormap(hot);
            caxis([min(min(ul1(:,:,kk))) max(max(ul1(:,:,kk)))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            xlabel('y [m]')
            ylabel('x [m]');
            title('u_l [m/s] Truncated');
            hold off;
            drawnow;
            
            if Movie==1; F = getframe(hfig); writeVideo(aviobj,F); end
        end
    end
    
    % Some error analysis between original and truncated sys
    RMSE(kk)               = rms(x(:,kk)-x1(:,kk));
    [maxe(kk),maxeloc(kk)] = max(abs(x(:,kk)-x1(:,kk)));
    
end

if Movie==1;close(aviobj);end
%%
figure(2);clf
subplot(2,1,1)

if Control>0
    if Truncated>0
        plot(time(1:end-1),znc1(1:end-1));hold on;
        plot(time(1:end-1),z1(1:end-1),'r');grid;
        title('\Sigma \delta u_r 2nd turbine \color{blue}uncontrolled, \color{red}controlled')
        ylabel('\Sigma u_r [m/s]');
    else
        plot(time(1:end-1),znc(1:end-1));hold on;
        plot(time(1:end-1),z(1:end-1),'r');grid;
        title('\Sigma \delta u_r 2nd turbine \color{blue}uncontrolled, \color{red}controlled')
        ylabel('\Sigma u_r [m/s]');
    end
else
    plot(time(1:end-1),znc(1:end-1));hold on;
    plot(time(1:end-1),znc1(1:end-1),'r');grid;
    title('\Sigma \delta u_r 2nd turbine \color{blue}Full, \color{red}Truncated')
    ylabel('\Sigma u_r [m/s]');
end

subplot(2,1,2)
plot(time(1:end-1),w(1,1:end-1));hold on;
plot(time(1:end-1),uc(1:end-1),'r');grid;
ylim([min(min(w(1,1:end-1)),min(uc(1:end-1)))-.1 max(max(w(1,1:end-1)),max(uc(1:end-1)))+.05]);
ylabel('\delta \beta');xlabel('Time [s]');
title('\color{blue}\delta \beta_1, \color{red}\delta \beta_2')