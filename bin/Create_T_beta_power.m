function controller = Create_T_beta_power(sys,Wp)

% This linear model uses the deviation power as output and beta as an
% input.

Em = sys.derivatives.A;
Am = blkdiag(spdiags(sys.derivatives.ccx,0,length(sys.derivatives.ccx),length(sys.derivatives.ccx)),...
    spdiags(sys.derivatives.ccy,0,length(sys.derivatives.ccy),length(sys.derivatives.ccy)),...
    0.*speye((Wp.mesh.Nx-2)*(Wp.mesh.Ny-2)));
Am = Am+sys.derivatives.dSmdx-sys.derivatives.dAdx;
Bm = sys.derivatives.dSmdbeta;
Cm = -sys.derivatives.dJdx';
Dm = -sys.derivatives.dJdbeta;

qq    = length(sys.derivatives.ccx)+length(sys.derivatives.ccy);

controller.model = dss(sys.Qsp'*full(Am(1:qq,1:qq))*sys.Qsp,sys.Qsp'*full(Bm(1:qq,:)),...
    full(Cm(:,1:qq))*sys.Qsp,Dm,sys.Qsp'*full(Em(1:qq,1:qq))*sys.Qsp,Wp.sim.h);

controller.A  = (sys.Qsp'*full(Em(1:qq,1:qq))*sys.Qsp)\(sys.Qsp'*full(Am(1:qq,1:qq))*sys.Qsp);
controller.B  = (sys.Qsp'*full(Em(1:qq,1:qq))*sys.Qsp)\(sys.Qsp'*full(Bm(1:qq,:)));
controller.C  = full(Cm(:,1:qq))*sys.Qsp;
controller.Cz = 1;
controller.D  = Dm;

% The model without direct feedthrough. 
controller.An = [controller.A controller.B; zeros(2,size(controller.A,1)+size(controller.B,2))];
controller.Bn = [zeros(size(controller.A,2),size(controller.B,2));eye(size(controller.B,2))];
controller.Cn = [controller.C controller.D];
controller.Dn = zeros(size(controller.Cn,1),size(controller.Bn,2));

controller.A  = controller.An;
controller.B  = controller.Bn;
controller.C  = controller.Cn;
controller.D  = controller.Dn;

%%
% model1=ss(full(controller.A),full(controller.B),full(controller.C),full(controller.D),1);
% model2=ss(full(controller.An),full(controller.Bn),full(controller.Cn),full(controller.Dn),1);
% 
% t  = (0:1:400)';
% r1 = 1*ones(length(t),1);
% r2 = 1*ones(length(t),1);
% r  = [r1 r2];
% 
% y1 = lsim(model1,r,t,zeros(size(model1.a,1),1));
% %y2 = lsim(model2,r,t,[zeros(size(model2.a,1)-2,1);Dm(1);Dm(2)]);
% y2 = lsim(model2,r,t,zeros(size(model2.a,1),1));
% 
% figure(1);clf;
% plot(t,y1);hold on;
% plot(t,y2);grid;