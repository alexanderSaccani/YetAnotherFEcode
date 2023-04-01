clearvars
close all
clc

%% Thermoelastic clamped plate

%Material _________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
alpha_T = 11.7e-6; % thermal expansion coefficient

myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'THERMAL_EXPANSION_COEFFICIENT', alpha_T);

%Geometry__________________________________________________________________
Lx = 0.4;
Ly = .25;
thickness = .001;     % thickness of plate
height_midspan = 4e-3; %height of midspan 

%Mesh______________________________________________________________________
nx = 15; % # elements along x
ny = 8; % # elements along y

eltype = 'QUAD8'; %'QUAD4'
myElementConstructor = @()Quad8Shell(thickness, myMaterial,2);
[nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,eltype);

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

nNodes = myMesh.nNodes;

%define z to assign to nodes
w = height_midspan/Lx;
R = (Lx/4 + Lx*w^2)/(2*w);
yy = R - Lx*w;
xx = Lx/2;
th0 = atan(yy/xx);

x = nodes(:,1);
zNodes = -R*sin(th0)+sqrt(-(x-Lx/2).^2+R^2);

%assign height z to nodes
for ii = 1:myMesh.nElements 
   myMesh.Elements(ii).Object.setZ(zNodes);   
end

%plot mesh
PlotMesh([nodes,zNodes],elements);

%Boundary conditions: all sides are clamped________________________________
myMesh.set_essential_boundary_condition([nset{1}],1:3,0)
myMesh.set_essential_boundary_condition([nset{2}],1:3,0)
myMesh.set_essential_boundary_condition([nset{3}],1:3,0)
myMesh.set_essential_boundary_condition([nset{4}],1:3,0)

%Assembly__________________________________________________________________
Assembly = Assembly(myMesh);
M = Assembly.mass_matrix();
u0 = zeros( myMesh.nDOFs, 1);

%Compute and store Mass matrix, K cold_____________________________________
%T field (zeros)
T = zeros(nNodes,1);
gradT = zeros(nNodes,1);
[Kc,~] = Assembly.tangent_stiffness_and_force(u0,T,gradT);

alpha = 0.01;
beta = 0.01;
D = alpha*Kc + beta*M;

% store matrices
Assembly.DATA.Kc = Kc; %Kcold
Assembly.DATA.M = M; %mass matrix
Assembly.DATA.C = D; %damping matrix


%% Vibration modes (cold)
% nVMs = 5; % first n_VMs modes with lowest frequency calculated 
% 
% Kconstr = Assembly.constrain_matrix(Assembly.DATA.Kc);
% Mconstr = Assembly.constrain_matrix(Assembly.DATA.M);
% 
% [VMsCold,om] = eigs(Kconstr, Mconstr, nVMs, 'SM');
% [fcold,ind] = sort(sqrt(diag(om))/2/pi);
% 
% VMsCold = VMsCold(:,ind);
% for ii = 1:nVMs
%     VMsCold(:,ii) = VMsCold(:,ii)/max(sqrt(sum(VMsCold(:,ii).^2,2)));
% end
% VMsCold = Assembly.unconstrain_vector(VMsCold);
% 
% %plot VMs__________________
% VMs = cell(nVMs,1);
% for ii = 1:nVMs  
%     VMs{ii,1} = reshape(VMsCold(:,ii),5,[]).';    
% end
% 
% for ii = 1:nVMs
%     figure; hold on;
%     VM2plot = VMs{ii};
%     PlotFieldonDeformedMesh([nodes,zNodes],elements, [VM2plot(:,1),...
%         VM2plot(:,2), VM2plot(:,3)], 'factor',max(nodes(:,2)));
%     title(['VM ',num2str(ii),', f: ', num2str(fcold(ii)),' Hz']);
% end

%% VMs and MDs Cold

T_ampls = 0;

n_sampl = length(T_ampls);
T_samples = zeros(myMesh.nNodes,n_sampl);
gradT_samples = zeros(myMesh.nNodes,n_sampl);

for ii = 1:n_sampl
    T_samples(:,ii) = ones(myMesh.nNodes,1)*T_ampls(ii);
end

nVMs = 4;

logic_MD = 1;
VMMD = VMMD_thermal(Assembly,nVMs,logic_MD,T_samples,gradT_samples);
fcold = VMMD.omega/2/pi;
VMsCold = VMMD.VMs;
MDsCold = VMMD.MDs.data;
MDsNames = VMMD.MDs.names;
eqPoint = VMMD.eq;

%Plot VMs
for jj = 1:length(T_ampls)
    
    VMs = cell(nVMs,1);
    for ii = 1:nVMs  
        VMs{ii,1} = reshape(VMsCold{jj}(:,ii),5,[]).';    
    end

    for ii = 1:nVMs
        figure; hold on;
        VM2plot = VMs{ii};
        PlotFieldonDeformedMesh([nodes,zNodes],elements, [VM2plot(:,1),...
            VM2plot(:,2), VM2plot(:,3)], 'factor', max(nodes(:,2)));
        title(['VM ',num2str(ii),' dT = ',num2str(T_ampls(jj)),'K, f: ', num2str(fcold(ii,jj)),' Hz']);
    end

end


% %plot MDs
% nMDs = size(MDsCold{1,1},2);
% 
% for jj = 1:length(T_ampls)
%     
%     MDs = cell(nMDs,1);
%     for ii = 1:nMDs  
%        MDs{ii,1} = reshape(MDsCold{jj}(:,ii),5,[]).';    
%     end
% 
%     for ii = 1:nMDs
%         figure; hold on;
%         MD2plot = MDs{ii};
%         PlotFieldonDeformedMesh([nodes,zNodes],elements, [MD2plot(:,1),...
%             MD2plot(:,2), MD2plot(:,3)], 'factor',max(nodes(:,2)));
%         title(['MD ',num2str(MDsNames(ii,1)),num2str(MDsNames(ii,2)),'  dT = ',num2str(T_ampls(jj)),'K']);
%     end
% 
% end


%% Buckling Analysis for different temperatures
% 
% BUCKLING_ANALYSIS = 1; %set to 0 to use data stored in variable "filename"
% filename = "bucklAnalFlat.mat";
% 
% if BUCKLING_ANALYSIS == 1
%     
%     %define forcing
%     F = zeros(myMesh.nDOFs,1);
%     nf = find_node(Lx/2,Ly/2,[],nodes); % node where to put the force
%     node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
%     F(node_force_dofs(3)) = +1;
% 
%     F0 = F; %load configuration corresponding to lambda = 1 (load multiplier)
% 
%     lam_high = 50; %load multiplier max abs value (symmetric analysis for positive and negative values of lambda from the origin
% 
%     %parameters for continuation
%     Sopt.parametrization = 'orthogonal';
%     Sopt.predictor = 'tangent';
%     Sopt.noconv_stepcor = 'red';
%     Sopt.errmax = 6;
%     Sopt.reversaltolerance = 1;
%     ds = 0.001;
%     Sopt.dsmin = ds/100;
%     Sopt.dsmax = ds*1000;
% 
%     %T samples 
%     T_samples = [0];
%     gradT = zeros(myMesh.nNodes,1);
% 
%     %initialize output
%     X_buckl_an = cell(length(T_samples),1);
% 
%     Solopt = optimset(optimset(@fsolve),'Display','off',...%'iter',...
%             'Jacobian','on','MaxIter',50,'DerivativeCheck','off');
% 
%     for ii = 1:length(T_samples)
% 
%         T_ampl_ii = T_samples(ii);
% 
%         T = T_ampl_ii*ones(myMesh.nNodes,1); %static temp profile
% 
%         [K_hot,gth] = Assembly.tangent_stiffness_and_force(zeros(myMesh.nDOFs,1),T,gradT);
%     %   u_lin_0 = Assembly.solve_system(K_hot, zeros(myMesh.nDOFs,1), gth);   %subtract the internal force generated  by the temperature
% 
%         u_lin_FAKE = Assembly.solve_system(Assembly.DATA.Kc, zeros(myMesh.nDOFs,1), gth); %compute the guess with cold tangent stiffness (otherwise it explodes)
%     %   u_guess = Assembly.constrain_vector(u_lin_0);
% 
%         u_guess = Assembly.constrain_vector(u_lin_FAKE);
%         % ds = norm(u_lin_0);
% 
%         fun_residual = @(X) residual_buckling(X,Assembly,F0,T,gradT);
%         [X1,Solinfo,Sol] = solve_and_continue(u_guess,fun_residual,0,lam_high,ds,Sopt,[],Solopt);
%         [X2,Solinfo,Sol] = solve_and_continue(u_guess,fun_residual,0,-lam_high,ds,Sopt,[],Solopt);
%         
%         X_buckl_an_ii = [flip(X2,2),X1];
%         lam_ii = X_buckl_an_ii(end,:);
%         X_ii = Assembly.unconstrain_vector(X_buckl_an_ii(1:end-1,:));
% 
%         X_buckl_an{ii} = [X_ii;lam_ii];
% 
% 
%     end
% 
% %      save("bucklAnal.mat","T_samples","X_buckl_an","F0");
% 
% else
%     
%     load(filename)
% 
% end
% 
% % plot results of buckling analysis
% plot_dof_location = [0.5,0.5]; %percentage of length of beam
% node2plot = find_node(plot_dof_location(1)*Lx,plot_dof_location(2)*Ly,[],nodes); % node to plot
% nDofPerNode = 5;
% dof2plot = get_index(node2plot, nDofPerNode);
% dof2plot = dof2plot(3); %plot vertical displacement
% 
% color_list = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE',...
%     '#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
% 
% figure('units','normalized','position',[.3 .3 .4 .4]);
% title(['buckling analysis: tranversal displacement at [', ...
%     num2str(plot_dof_location(1)*100),',',num2str(plot_dof_location(2)*100),']','% of panel span']);
% leg = cell(length(T_samples),1); %legend initialize
% for ii = 1:length(T_samples)
%     
%     hold on;
%     X_ii = X_buckl_an{ii};
%     delta = X_ii(dof2plot,:);
%     lam = X_ii(end,:);
%     plot(delta,lam,'color',color_list{ii},'linewidth',2);
%     
%     xlabel('displacement [m]')
%     ylabel('\lambda')
%     
%     leg{ii} = ['T = ',num2str(T_samples(ii))];
%     
% end
% legend(leg);


%% Hot modes

% T_ampls = [5,10];
% 
% n_sampl = length(T_ampls);
% T_samples = zeros(myMesh.nNodes,n_sampl);
% gradT_samples = zeros(myMesh.nNodes,n_sampl);
% 
% for ii = 1:n_sampl
%     T_samples(:,ii) = ones(myMesh.nNodes,1)*T_ampls(ii);
% end
% 
% nVMs = 5;
% [VMsHot,static_eq,omega] = VM_thermal(Assembly,nVMs,T_samples,gradT_samples);
% fhot = omega/2/pi;
% 
% %Plots
% for jj = 1:length(T_ampls)
%     
%     VMs = cell(nVMs,1);
%     for ii = 1:nVMs  
%         VMs{ii,1} = reshape(VMsHot{jj}(:,ii),5,[]).';    
%     end
% 
%     for ii = 1:nVMs
%         figure; hold on;
%         VM2plot = VMs{ii};
%         PlotFieldonDeformedMesh([nodes,zNodes],elements, [VM2plot(:,1),...
%             VM2plot(:,2), VM2plot(:,3)], 'factor',max(nodes(:,2)));
%         title(['VM ',num2str(ii),' dT = ',num2str(T_ampls(jj)),'K, f: ', num2str(fhot(ii,jj)),' Hz']);
%     end
% 
% end

%% Compute Hot modes and Modal Derivatives

T_ampls = [5,10,15];

n_sampl = length(T_ampls);
T_samples = zeros(myMesh.nNodes,n_sampl);
gradT_samples = zeros(myMesh.nNodes,n_sampl);

for ii = 1:n_sampl
    T_samples(:,ii) = ones(myMesh.nNodes,1)*T_ampls(ii);
end

nVMs = 4;

logic_MD = 1;
VMMD = VMMD_thermal(Assembly,nVMs,logic_MD,T_samples,gradT_samples);
fhot = VMMD.omega/2/pi;
VMsHot = VMMD.VMs;
MDsHot = VMMD.MDs.data;
MDsNames = VMMD.MDs.names;
eqPoint = VMMD.eq;

% %Plot equilibrium configuration
% for jj = 1:length(T_ampls)
%     eqPoint = reshape(eqPoint,5,[]).';
%     figure
%     PlotFieldonDeformedMesh([nodes,zNodes],elements, [eqPoint(:,1),...
%             eqPoint(:,2), eqPoint(:,3)], 'factor',5*max(nodes(:,2)));
%         title(['equilibrium position for uniform dT = ',num2str(T_ampls(jj)),'K']);
% end
% 
%Plot VMs
for jj = 1:length(T_ampls)
    
    VMs = cell(nVMs,1);
    for ii = 1:nVMs  
        VMs{ii,1} = reshape(VMsHot{jj}(:,ii),5,[]).';    
    end

    for ii = 1:nVMs
        figure; hold on;
        VM2plot = VMs{ii};
        PlotFieldonDeformedMesh([nodes,zNodes],elements, [VM2plot(:,1),...
            VM2plot(:,2), VM2plot(:,3)], 'factor',5*max(nodes(:,2)));
        title(['VM ',num2str(ii),' dT = ',num2str(T_ampls(jj)),'K, f: ', num2str(fhot(ii,jj)),' Hz']);
    end

end
% 
% 
% %plot MDs
% nMDs = size(MDsHot{1,1},2);
% 
% for jj = 1:length(T_ampls)
%     
%     MDs = cell(nMDs,1);
%     for ii = 1:nMDs  
%        MDs{ii,1} = reshape(MDsHot{jj}(:,ii),5,[]).';    
%     end
% 
%     for ii = 1:nMDs
%         figure; hold on;
%         MD2plot = MDs{ii};
%         PlotFieldonDeformedMesh([nodes,zNodes],elements, [MD2plot(:,1),...
%             MD2plot(:,2), MD2plot(:,3)], 'factor',0.1*max(nodes(:,2)));
%         title(['MD ',num2str(MDsNames(ii,1)),num2str(MDsNames(ii,2)),'  dT = ',num2str(T_ampls(jj)),'K']);
%     end
% 
% end



%% Define dynamic temperature variation

T0 = 1*ones(myMesh.nNodes,1);
Tinitial = 0;
Tfinal = 12;
eps = 0.1; 
timeRamp = 1/fcold(1)/eps;

Tdynt = @(t) T0*(Tfinal - Tinitial)/timeRamp*t + (T0*Tfinal - T0*(Tfinal - Tinitial)/timeRamp*t)*heaviside(t-timeRamp);
gradTt = @(t) zeros(myMesh.nNodes,1);

%% Define forcing
F0 = zeros(myMesh.nDOFs,1);

nf = find_node(Lx/3,Ly/3,[],myMesh.nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F0(node_force_dofs(3)) = +150;

fForc = (fcold(1)+fcold(2))/2;
omForc = 2*pi*fForc;
TForc = 1/fForc;

Fext = @(t) F0*cos(omForc*t);

%% Define integration time
tint = 3*timeRamp;

%% Integrate full model
% settings for integration
h = 2*pi/omForc/50; % time step for integration

% Initial condition: equilibrium
nDofFull = length(Assembly.Mesh.EBC.unconstrainedDOFs);

q0 = zeros(nDofFull,1);
qd0 = zeros(nDofFull,1);
qdd0 = zeros(nDofFull,1);

% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear Residual evaluation function handle
residualFull = @(q,qd,qdd,t)residual_thermal_nonlinear(q,qd,qdd,t,Assembly,Fext,Tdynt,gradTt);

% integrate equations with Newmark scheme
TI_NL.Integrate(q0,qd0,qdd0,tint,residualFull);

% project to full space
uDynFull = TI_NL.Solution.q;

%save solution
uDynFull = Assembly.unconstrain_vector(uDynFull); 
dynFull.nlin.disp = decodeDofsNodes(uDynFull,nNodes,5); % (node, dof of node, tsamp)
dynFull.nlin.time = TI_NL.Solution.time;


% plot_dof_location = [0.5,0.5]; %percentage of length of beam
% node2plot = find_node(plot_dof_location(1)*Lx,plot_dof_location(2)*Ly,[],nodes); % node to plot
% 
% figure
% plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,3,:)))

%% Construct ROM
%compute the static solution (to forcing)
[~,uStatic] = static_eq_shell_thermal(Assembly, F0, zeros(myMesh.nNodes,1), zeros(myMesh.nNodes,1), 'display', 'iter-detailed');
% UNL = reshape(u,5,[]).';        % Nonlinear response

Vcold = orth([VMsCold{1,1},MDsCold{1,1},uStatic]); %reduction basis (cold basis)

V = Vcold;

redAssembly = ReducedAssembly(myMesh,V); %reduced assembly
redAssembly.DATA.M = V.'*Assembly.DATA.M*V;
redAssembly.DATA.C = V.'*Assembly.DATA.C*V;

Tred = zeros(nNodes,1);
gradTred = zeros(nNodes,1);
%[Ktgred,Fred] = RedAssembly.tangent_stiffness_and_force(zeros(size(Vcold,2),1),Tred,gradTred);

%% Integrate nonlinear ROM (constant basis)
% settings for integration
h = 2*pi/omForc/50; % time step for integration

% Initial condition: equilibrium
nDofRed = size(V,2);

q0 = zeros(nDofRed,1);
qd0 = zeros(nDofRed,1);
qdd0 = zeros(nDofRed,1);

% Instantiate object for nonlinear time integration
TI_NL_r = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear Residual evaluation function handle
residual_r = @(q,qd,qdd,t)residual_thermal_reduced_nonlinear(q,qd,qdd,t,redAssembly,Fext,Tdynt,gradTt);

% integrate equations with Newmark scheme
TI_NL_r.Integrate(q0,qd0,qdd0,tint,residual_r);

% project to full space
uDynRed = Vcold*TI_NL_r.Solution.q;

%save solution
dynRed.nlin.disp = decodeDofsNodes(uDynRed,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dynRed.nlin.time = TI_NL_r.Solution.time;


%% plot comparison
plotDofLocation = [0.3,0.5]; %percentage of plate length (along x and y)
node2plot = find_node(plotDofLocation(1)*Lx,plotDofLocation(2)*Ly,[],nodes); % node to plot

fontsize = 15;
linewidth = 2;
figure('units','normalized','position',[0.3,0.3,0.4,0.4]); hold on
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,3,:)),'-k','linewidth',linewidth)
plot(dynRed.nlin.time,squeeze(dynRed.nlin.disp(node2plot,3,:)),'--r','linewidth',linewidth)
xlabel('time [s]','fontsize',fontsize); ylabel('u [m]','fontsize', fontsize)
ax = gca; grid on; ax.FontSize = fontsize; legend('HFM','RED','fontsize',fontsize)

