function controller = Computecontrolsignal(Wp,sys,controller,sol,perturbatie,scriptOptions)

global sysMpc model MV mpcobj xmpc mpcoptions uc

beta       = Wp.turbine.input(sol.k).beta;
phi        = Wp.turbine.input(sol.k).phi;
CT_prime   = Wp.turbine.input(sol.k).CT_prime;
Rho        = Wp.site.Rho;
Drotor     = Wp.turbine.Drotor;
powerscale = Wp.turbine.powerscale;
Ar         = pi*(0.5*Drotor)^2;
constant   = powerscale*.5*Rho*Ar;

if sol.k>=1
    
    % linearisation points
    controller.ss.U         = zeros(Wp.turbine.N,1);
    uss                     = sol.u;
    for kk=1:Wp.turbine.N
        controller.ss.U(kk) = sum(uss(Wp.mesh.xline(kk),Wp.mesh.yline{kk}))/length(Wp.mesh.yline{kk});
    end
    controller.ss.P         = sol.turbine.power;
    controller.ss.CT_prime  = sol.turbine.CT_prime;
    
    % System matrices
    if scriptOptions.Projection>0
        controller.A          = sys.Etl\sys.Atl;
        nx                    = size(controller.A,1);
        controller.B          = sys.Etl\sys.Btl;
        controller.Cz         = zeros(Wp.turbine.N-1,size(sys.Qsp,1));
        for kk=1:Wp.turbine.N
            controller.Cz(kk,(Wp.mesh.xline(kk)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{kk}(1)-1:...
                (Wp.mesh.xline(kk)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{kk}(end)-1) = 1/length(Wp.mesh.yline{kk});
        end
        controller.C          = sys.Qsp;
        controller.D          = zeros(size(controller.Cz,1),size([beta;phi],1));
        ny                    = size(controller.C,1);
        
    else
        % sys.A * dx_{k+1} = sys.Al * dx_k + sys.Bl * du_k
        controller.A          = sys.A\sys.Al;
        nx                    = size(controller.A,1);
        controller.B          = sys.A\sys.Bl;
        nu                    = size(controller.B,2);
        controller.Cz         = zeros(Wp.turbine.N-1,size(sol.x,1));
        for kk=1:Wp.turbine.N
            controller.Cz(kk,(Wp.mesh.xline(kk)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{kk}(1)-1:...
                (Wp.mesh.xline(kk)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{kk}(end)-1) = 1/length(Wp.mesh.yline{kk});
        end
        controller.C          = sparse(eye(size(sol.x,1)));
        controller.Cx         = diag(3*constant*controller.ss.U.^2.*controller.ss.CT_prime);
        controller.Du         = [diag(constant*controller.ss.U.^3) sparse(Wp.turbine.N,Wp.turbine.N)];
        controller.D          = zeros(size(controller.Cz,1),size([CT_prime;phi],1));
        ny                    = size(controller.Cz,1);
    end
    
    if scriptOptions.Control
        sysMpc     = (ss(full(blkdiag(controller.A,0*eye(2*Wp.turbine.N))),[full(controller.B);eye(nu,nu)],...
            full([controller.Cx(2:end,:)*controller.Cz controller.Du(2:end,:)]),...
            controller.D(2:end,:),Wp.sim.h));
        
        if Wp.turbine.N==2
            model      = setmpcsignals(sysMpc,'MD',[1 3 4],'MV',2,'MO',1);
            MV         = struct('Min',-.25,'Max',.6,'RateMin',-.01,'RateMax',.01);
        elseif Wp.turbine.N==3
            model      = setmpcsignals(sysMpc,'MD',[1 4 5 6],'MV',[2 3],'MO',[1 2]);
            MV         = struct('Min',{-.45;-.45},'Max',{.5;.5},'RateMin',{-.1;-.1},'RateMax',{.1;.1});
        elseif Wp.turbine.N==6
            model      = setmpcsignals(sysMpc,'MD',[1 7 8 9 10 11 12],'MV',[2 3 4 5 6],'MO',[1 2 3 4 5]);
            MV         = struct('Min',{-1;-1;-1;-1;-1},'Max',{1;1;1;1;1},'RateMin',{-.1;-.1;-.1;-.1;-.1},'RateMax',{.1;.1;.1;.1;.1});
        end
        
        np = 250;       % Prediction horizon
        nc = 20;        % Control horizon
        
        mpcobj     = mpc(model,Wp.sim.h,np,nc,[],MV);
        xmpc       = mpcstate(mpcobj);
        mpcoptions = mpcmoveopt;
    end
    
    if sol.k==1
        % vectors for simulation
        controller.u   = zeros(2*Wp.turbine.N,Wp.sim.NN);
        controller.x   = zeros(nx,Wp.sim.NN);
        controller.y   = zeros(ny,Wp.sim.NN);
        controller.zl  = zeros(Wp.turbine.N,Wp.sim.NN);         % Averaged rotor velocities linear model
        controller.znl = zeros(Wp.turbine.N,Wp.sim.NN);
        controller.z   = zeros(Wp.turbine.N,Wp.sim.NN);         % Averaged rotor velocities nonlinear model
        
        % reference
        %controller.r   = repmat((Ur-Uo),1,Wp.sim.NN);          % Reference to sum of rotor velocity of downwind turbines
        controller.r   = zeros(Wp.turbine.N-1,Wp.sim.NN);       % For disturbance rejection
    end
    
end

%%
controller.znl(:,sol.k)  = sol.turbine.power;                                     % nonlinear model/measurement

% Compute control downwind turbines
if scriptOptions.Control
    uc                   = mpcmove(mpcobj,xmpc,controller.znl(2:end,sol.k)-controller.ss.P(2:end),controller.r(:,sol.k),...
        [perturbatie;phi],mpcoptions);
else
    uc                   = zeros(Wp.turbine.N-1,1);
end

controller.u(:,sol.k)    = [perturbatie;uc;phi];

controller.x(:,sol.k+1)  = controller.A*controller.x(:,sol.k) + controller.B*controller.u(:,sol.k);
controller.y(:,sol.k)    = controller.Cz*controller.x(:,sol.k);                         % deviation variable dU
if sol.k==1
    controller.z(:,sol.k)    = [controller.Cx controller.Du]*...
        [controller.y(:,sol.k);controller.u(:,sol.k)];      % deviation variable dP
else
    controller.z(:,sol.k)    = [controller.Cx controller.Du]*...
        [controller.y(:,sol.k);controller.u(:,sol.k-1)];    % deviation variable dP
end
controller.zl(:,sol.k)   = controller.z(:,sol.k) + controller.ss.P;                     % linear model

end