function controller = Computecontrolsignal(Wp,sys,controller,sol,input,k,perturbatie)

global sysMpc model MV mpcobj xmpc mpcoptions uc

if k==3
    
    % System matrices
    controller.A          = sys.Etl\sys.Atl;
    nx                    = size(controller.A,1);
    controller.B          = sys.Etl\sys.Btl; 
    controller.Cz         = zeros(Wp.turbine.N-1,size(sys.Qsp,1));
    for kk=2:Wp.turbine.N
        controller.Cz(kk-1,(Wp.mesh.xline(kk)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{kk}(1)-1:...
            (Wp.mesh.xline(kk)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{kk}(end)-1) = 1/length(Wp.mesh.yline{kk});
    end
    controller.C          = sys.Qsp;
    controller.D          = zeros(size(controller.Cz,1),size([input.beta;input.phi],1));
    ny                    = size(controller.C,1);
    
    sysMpc     = ss(full(controller.A),full(controller.B),full(controller.Cz*controller.C),controller.D,Wp.sim.h);
    
    if Wp.turbine.N==2
    model      = setmpcsignals(sysMpc,'MD',[1 3 4],'MV',2,'MO',1);
    MV         = struct('Min',-.25,'Max',.6,'RateMin',-.01,'RateMax',.01);
    elseif Wp.turbine.N==3
    model      = setmpcsignals(sysMpc,'MD',[1 4 5 6],'MV',[2 3],'MO',[1 2]);
    MV         = struct('Min',{-.45;-.45},'Max',{.5;.5},'RateMin',{-.05;-.05},'RateMax',{.05;.05});
    end 
    
    np = 100;   % Prediction horizon
    nc = 1;     % Control horizon
    
    mpcobj     = mpc(model,Wp.sim.h,np,nc,[],MV);
    
    if k==3
        controller.x       = zeros(nx,Wp.sim.NN);
        controller.y       = zeros(ny,Wp.sim.NN);
        controller.ss      = zeros(Wp.turbine.N-1,1);
        uss                = sol.u;
        
        for kk=2:Wp.turbine.N
            controller.ss(kk-1) = sum(uss(Wp.mesh.xline(kk),Wp.mesh.yline{kk}))/length(Wp.mesh.yline{kk});
        end
        
        controller.zl      = zeros(Wp.turbine.N-1,Wp.sim.NN);       % Normalised rotor velocities linear model
        controller.znl     = zeros(Wp.turbine.N-1,Wp.sim.NN);
        controller.z       = zeros(Wp.turbine.N-1,Wp.sim.NN);       % Normalised rotor velocities nonlinear model
        controller.r       = 0*ones(Wp.turbine.N-1,Wp.sim.NN);      % Reference to sum of rotor velocity of downwind turbines
    end
    
    
    if k==3 % xp,xd,xn,u,p
        xmpc      = mpcstate(mpcobj);
    end
    mpcoptions    = mpcmoveopt;
    
end

controller.y(:,k)    = controller.C*controller.x(:,k);
controller.z(:,k)    = controller.Cz*controller.y(:,k);                     % Sum rotor velocities turbine 2 deviation variable
controller.zl(:,k)   = controller.z(:,k) + controller.ss;
controller.znl(:,k)  = controller.Cz*sol.x;                                 % Sum rotor velocities turbine 2 nonlinear system

% Compute control downwind turbines
uc                   = mpcmove(mpcobj,xmpc,controller.znl(:,k)-controller.ss,controller.r(:,k),[perturbatie;input.phi],mpcoptions);
uc                   = zeros(Wp.turbine.N-1,1); %  Uncomment for open-loop
controller.x(:,k+1)  = controller.A*controller.x(:,k) + controller.B*[perturbatie;uc;input.phi];


end