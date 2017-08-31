function controller = Computecontrolsignal(Wp,sys,controller,sol,input,k,perturbatie,options)

global sysMpc model MV mpcobj xmpc mpcoptions uc

if k==3
    
    % System matrices
    if options.Projection>0
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
        
    else
        % sys.A * dx_{k+1} = sys.Al * dx_k + sys.Bl * du_k
        controller.A          = sys.A\sys.Al;
        nx                    = size(controller.A,1);
        controller.B          = sys.A\sys.Bl;
        controller.Cz         = zeros(Wp.turbine.N-1,size(sol.x,1));
        for kk=2:Wp.turbine.N
            controller.Cz(kk-1,(Wp.mesh.xline(kk)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{kk}(1)-1:...
                (Wp.mesh.xline(kk)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{kk}(end)-1) = 1/length(Wp.mesh.yline{kk});
        end
        controller.C          = sparse(eye(size(sol.x,1)));
        controller.D          = zeros(size(controller.Cz,1),size([input.beta;input.phi],1));
        ny                    = size(controller.C,1);
    end
    
    sysMpc     = ss(full(controller.A),full(controller.B),full(controller.Cz*controller.C),controller.D,Wp.sim.h);
    
    if Wp.turbine.N==2
        model      = setmpcsignals(sysMpc,'MD',[1 3 4],'MV',2,'MO',1);
        MV         = struct('Min',-.25,'Max',.6,'RateMin',-.01,'RateMax',.01);
    elseif Wp.turbine.N==3
        model      = setmpcsignals(sysMpc,'MD',[1 4 5 6],'MV',[2 3],'MO',[1 2]);
        MV         = struct('Min',{-.45;-.45},'Max',{.5;.5},'RateMin',{-.1;-.1},'RateMax',{.1;.1});
    end
    
    np = 250;   % Prediction horizon
    nc = 100;     % Control horizon
    
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
        % Give Power reference
        Uo  = controller.ss;        
        CTo = controller.CT(2:3);
        Po  = controller.Power(2:3); 
        
        Pr  = Po+[-35000;-35000];
        Ur  = (Pr./(.5*Wp.site.Rho*pi*(0.5*Wp.turbine.Drotor)^2*.5*CTo)).^(1/3);
        
        controller.Pr = repmat(Pr,1,Wp.sim.NN);
        %controller.r  = repmat((Ur-Uo),1,Wp.sim.NN);  % Reference to sum of rotor velocity of downwind turbines
        controller.r  = 0*ones(Wp.turbine.N-1,Wp.sim.NN); % For disturbance rejection     
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

% uncomment for no control
%uc                   = zeros(Wp.turbine.N-1,1);

% uncomment for steady-state control (only two-turbine case)
% if k<=3
%     uc                   = zeros(Wp.turbine.N-1,1);
% else
%     uc                   = [0.3902;-0.0740];
% end

controller.x(:,k+1)  = controller.A*controller.x(:,k) + controller.B*[perturbatie;uc;input.phi];


end