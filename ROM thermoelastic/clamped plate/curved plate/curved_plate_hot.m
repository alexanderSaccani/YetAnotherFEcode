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

nDofPerNode = 5;

elementsPlot = elements(:,1:4); %for plots

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
figure
PlotMesh([nodes,zNodes],elementsPlot);
axis on

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
%     PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [VM2plot(:,1),...
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

nVMs = 5;

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
        PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [VM2plot(:,1),...
            VM2plot(:,2), VM2plot(:,3)], 'factor', 15*max(nodes(:,2)));
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
%         PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [MD2plot(:,1),...
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
%         PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [VM2plot(:,1),...
%             VM2plot(:,2), VM2plot(:,3)], 'factor',max(nodes(:,2)));
%         title(['VM ',num2str(ii),' dT = ',num2str(T_ampls(jj)),'K, f: ', num2str(fhot(ii,jj)),' Hz']);
%     end
% 
% end


%% Define dynamic temperature variation

%define the shape of T distribution

%1)constant T
% Tshape = 1*ones(myMesh.nNodes,1);

%2)Gaussian temperature
%independent parameters of T distribution
sigma3X = Lx/2;
sigma3Y = Ly/2;

theta = 0; %inclination w.r.t. x axis
xc = Lx/3; %center of pulse
yc = Ly/3; %center of pulse

%dependent parameteres
sigmaX = sigma3X/3;
sigmaY = sigma3Y/3;

a = cos(theta)^2/(2*sigmaX^2) + sin(theta)^2/(2*sigmaY^2);
b = - sin(2*theta)/(4*sigmaX^2) + sin(2*theta)/(4*sigmaY^2);
c = sin(theta)^2/(2*sigmaX^2) + cos(theta)^2/(2*sigmaY^2);

nodesX = nodes(:,1);
nodesY = nodes(:,2);

Tshape = exp( -( a*(nodesX-xc).^2 +2*b*(nodesX-xc).*(nodesY-yc) +c*(nodesY-yc).^2) ); %gaussian with center at xc,yc inclined of angle theta w.r.t. x axis


%plot Tshape
figure
hf = PlotFieldonMesh(nodes,elementsPlot,Tshape);
title('Temperature distribution shape')

Tinitial = 0;
Tfinal = 70;
eps = 0.1; 
timeRamp = 1/fcold(1)/eps;

tstartRamp = timeRamp;
%Adynt = @(t) (Tfinal - Tinitial)/timeRamp*t + (Tfinal - (Tfinal - Tinitial)/timeRamp*t).*heaviside(t-timeRamp);
Adynt = @(t) (Tfinal - Tinitial)/timeRamp*t + (Tfinal - (Tfinal - Tinitial)/timeRamp*t).*heaviside(t-timeRamp) + ...
         -heaviside(-t).*(Tfinal - Tinitial)/timeRamp.*t;
Adynt = @(t) Adynt(t-tstartRamp);

Tdynt = @(t) Tshape*(Adynt(t)).';
%Tdynt = @(t) (Tshape*(Tfinal - Tinitial)/timeRamp*(t-tstartRamp) + (Tshape*Tfinal - Tshape*(Tfinal - Tinitial)/timeRamp*(t-tstartRamp))*heaviside((t-tstartRamp)-timeRamp));

gradTt = @(t) zeros(myMesh.nNodes,1);

%% Compute Hot modes and Modal Derivatives
T_ampls = 70;

n_sampl = length(T_ampls);
T_samples = zeros(myMesh.nNodes,n_sampl);
gradT_samples = zeros(myMesh.nNodes,n_sampl);

for ii = 1:n_sampl
    T_samples(:,ii) = Tshape*T_ampls(ii);
end

nVMs = 5;

logic_MD = 1;
VMMD = VMMD_thermal(Assembly,nVMs,logic_MD,T_samples,gradT_samples);
fhot = VMMD.omega/2/pi;
VMsHot = VMMD.VMs;
MDsHot = VMMD.MDs.data;
MDsNames = VMMD.MDs.names;
eqPoint = VMMD.eq;

%Plot equilibrium configuration
for jj = 1:length(T_ampls)
    eqPoint = reshape(eqPoint,5,[]).';
    figure
    PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [eqPoint(:,1),...
            eqPoint(:,2), eqPoint(:,3)], 'factor',120*max(nodes(:,2)));
        %title(['equilibrium position for uniform dT = ',num2str(T_ampls(jj)),'K']);
end

%Plot VMs
for jj = 1:length(T_ampls)
    
    VMs = cell(nVMs,1);
    for ii = 1:nVMs  
        VMs{ii,1} = reshape(VMsHot{jj}(:,ii),5,[]).';    
    end

    for ii = 1:nVMs
        figure; hold on;
        VM2plot = VMs{ii};
        PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [VM2plot(:,1),...
            VM2plot(:,2), VM2plot(:,3)], 'factor',15*max(nodes(:,2)));
        title(['VM ',num2str(ii),' dT = ',num2str(T_ampls(jj)),'K, f: ', num2str(fhot(ii,jj)),' Hz']);
        colorbar off
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
%         PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [MD2plot(:,1),...
%             MD2plot(:,2), MD2plot(:,3)], 'factor',0.1*max(nodes(:,2)));
%         title(['MD ',num2str(MDsNames(ii,1)),num2str(MDsNames(ii,2)),'  dT = ',num2str(T_ampls(jj)),'K']);
%     end
% 
% end
%% Define forcing
% F0 = zeros(myMesh.nDOFs,1);

% nf = find_node(Lx/3,Ly/3,[],myMesh.nodes); % node where to put the force
% node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
% F0(node_force_dofs(3)) = +260;

%pressure (triangular)
nodesX = nodes(:,1);
Fmax = 5000/length(nodes);
Fmin = 0;
Fz = Fmax-nodesX*Fmax/Lx;
F0 = zeros(myMesh.nDOFs,1);
F0(3:5:end) = Fz;

%plot forcing shape
figure
hf = PlotFieldonMesh(nodes,elementsPlot,Fz);
title('Forcing distribution shape (out of plane)')

fForc = 160; %tra la seconda e la terza calde, tra la prima e la seconda fredde
omForc = 2*pi*fForc;
TForc = 1/fForc;

Fext = @(t) F0*cos(omForc*t);

%% Define integration time
tint = 3.5*timeRamp;

%% Plot dynamic T over the analysis
tPlot = linspace(0,tint,100);
figure
for ii = 1:length(tPlot)
    tii = tPlot(ii);
    hf = PlotFieldonMesh(nodes,elementsPlot,Tdynt(tii));
    title('Temperature distribution in time')
    clim([0, Tfinal]); 
    txt = text('Units', 'Normalized', 'Position', [0.5, 0.9], 'string',...
    ['t: ',num2str(round(tii,5)),'s'], 'fontsize', 12); drawnow limitrate
    delete(txt)
end
txt = text('Units', 'Normalized', 'Position', [0.4, 0.9], 'string',...
    ['t: ',num2str(round(tii,5)),'s'], 'fontsize', 12);

tt = linspace(0,tint,1000);
figure
plot(tt,Adynt(tt),'-k','linewidth',2)
xlabel('time [s]')
ylabel('T_{max} [K]')
%% Integrate full model
% settings for integration
h = 2*pi/omForc/25; % time step for integration

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

%% integrate linearized model
% settings for integration
h = 2*pi/omForc/25; % time step for integration

% Initial condition: equilibrium
nDofFull = length(Assembly.Mesh.EBC.unconstrainedDOFs);

q0 = zeros(nDofFull,1);
qd0 = zeros(nDofFull,1);
qdd0 = zeros(nDofFull,1);

% Instantiate object for nonlinear time integration
TI_LIN = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% nonlinear Residual evaluation function handle
residualFullLin = @(q,qd,qdd,t)residual_thermal_linear(q,qd,qdd,t,Assembly,Fext,Tdynt,gradTt);

% integrate equations with Newmark scheme
TI_LIN.Integrate(q0,qd0,qdd0,tint,residualFullLin);

% project to full space
uDynFullLin = TI_LIN.Solution.q;

%save solution
uDynFullLin = Assembly.unconstrain_vector(uDynFullLin); 
dynFull.lin.disp = decodeDofsNodes(uDynFullLin,nNodes,5); % (node, dof of node, tsamp)
dynFull.lin.time = TI_LIN.Solution.time;

plot_dof_location = [0.3,0.3]; %percentage of length of beam
node2plot = find_node(plot_dof_location(1)*Lx,plot_dof_location(2)*Ly,[],nodes); % node to plot

figure
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,3,:)))
hold on
plot(dynFull.lin.time,squeeze(dynFull.lin.disp(node2plot,3,:)))
legend('NL','LIN')

%% Construct ROM
%compute the static solution (to forcing)
[~,uStaticF] = static_eq_shell_thermal(Assembly, F0, zeros(myMesh.nNodes,1), zeros(myMesh.nNodes,1), 'display', 'iter-detailed');

%compute the static solution (to temperature)
argTanStiff = {Tfinal*Tshape, zeros(myMesh.nNodes,1)};
[~,Fth] = Assembly.tangent_stiffness_and_force(0*F0,argTanStiff{:});
Assembly.DATA.K = Assembly.DATA.Kc; %cold T used to predict the response with linear analysis
Assembly.DATA.F0 = Fth;
[~,uStaticT] = static_eq(Assembly, 0*F0, 'vararginTanStiffForce',argTanStiff, 'display', 'iter-detailed');

%compute the static solution (to forcing at hot T)
Assembly.DATA.K = Assembly.DATA.Kc; %cold T used to predict the response with linear analysis
Assembly.DATA.F0 = Fth;
[~,uStaticTF] = static_eq(Assembly, F0,  'vararginTanStiffForce',argTanStiff, 'display', 'iter-detailed');

% UNL = reshape(u,5,[]).';        % Nonlinear response

%Vcold = orth([VMsCold{1,1},MDsCold{1,1},uStaticF,uStaticT]);%,uStaticT]);%,uStaticTF]); %reduction basis (cold basis)
Vcold = orth([VMsCold{1,1},MDsCold{1,1}]);
Vcold = orth([VMsCold{1,1},MDsCold{1,1},uStaticF,uStaticT,uStaticTF]);
Vhot = orth([VMsHot{1,1},MDsHot{1,1},uStaticF,uStaticT,uStaticTF]);
%Vcold = orth([VMsCold{1,1}]);

V = Vcold;
%V = Vhot;

redAssembly = ReducedAssembly(myMesh,V); %reduced assembly
redAssembly.DATA.M = V.'*Assembly.DATA.M*V;
redAssembly.DATA.C = V.'*Assembly.DATA.C*V;

% Tred = zeros(nNodes,1);
% gradTred = zeros(nNodes,1);
%[Ktgred,Fred] = RedAssembly.tangent_stiffness_and_force(zeros(size(Vcold,2),1),Tred,gradTred);

%% Integrate nonlinear ROM (constant basis)
% settings for integration
h = 2*pi/omForc/25; % time step for integration

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
uDynRed = V*TI_NL_r.Solution.q;
vDynRed = V*TI_NL_r.Solution.qd; 
aDynRed = V*TI_NL_r.Solution.qdd;

%save solution
dynRed.nlin.disp = decodeDofsNodes(uDynRed,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dynRed.nlin.time = TI_NL_r.Solution.time;

%% compute error
nodesErr = 1:1:size(nodes,1);
gdlErr = 1:5;

%RMS error
errROM = RMS_error_ROM(dynFull.nlin.disp,dynFull.nlin.time,dynRed.nlin.disp,dynRed.nlin.time,gdlErr,nodesErr);
errROMrotation = RMS_error_ROM(dynFull.nlin.disp,dynFull.nlin.time,dynRed.nlin.disp,dynRed.nlin.time,4:5,nodesErr);
errROMwDisp = RMS_error_ROM(dynFull.nlin.disp,dynFull.nlin.time,dynRed.nlin.disp,dynRed.nlin.time,3,nodesErr);
errROMmembrane = RMS_error_ROM(dynFull.nlin.disp,dynFull.nlin.time,dynRed.nlin.disp,dynRed.nlin.time,1:2,nodesErr);

%err Res
ROMdisp = Assembly.constrain_vector(uDynRed);
ROMvel = Assembly.constrain_vector(vDynRed);
ROMacc = Assembly.constrain_vector(aDynRed);
ROMtime = dynRed.nlin.time;
errRes = residual_error_ROM(residualFull,ROMdisp,ROMvel,ROMacc,ROMtime);

%% plot comparison (displacements)
%displacements
plotDofLocation = [0.3,0.3]; %percentage of plate length (along x and y)
node2plot = find_node(plotDofLocation(1)*Lx,plotDofLocation(2)*Ly,[],nodes); % node to plot
dofPl = 3;

description = 'reduction basis: first 5 cold VMs and MDs + static T + static TF, nDofRed = 22 ';
fontsize = 15;
linewidth = 2;

figure('units','normalized','position',[0.3,0.1,0.6,0.7]); hold on
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,dofPl,:)),'-k','linewidth',linewidth)
plot(dynRed.nlin.time,squeeze(dynRed.nlin.disp(node2plot,dofPl,:)),'--r','linewidth',linewidth)
plot(dynFull.lin.time,squeeze(dynFull.lin.disp(node2plot,dofPl,:)),'--b','linewidth',linewidth)
xlabel('time [s]','fontsize',fontsize); ylabel('u [m]','fontsize', fontsize)
title(['normalized plotted dof location: ',num2str(plotDofLocation(1)),',',...
    num2str(plotDofLocation(2)),',  dof: ',num2str(dofPl)]);
ax = gca; grid on; ax.FontSize = fontsize; legend('HFM','RED','LIN','fontsize',fontsize);


figure('units','normalized','position',[0.3,0.1,0.6,0.7]); 
subplot 311; hold on
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,dofPl,:)),'-k','linewidth',linewidth)
plot(dynRed.nlin.time,squeeze(dynRed.nlin.disp(node2plot,dofPl,:)),'--r','linewidth',linewidth)
plot(dynFull.lin.time,squeeze(dynFull.lin.disp(node2plot,dofPl,:)),'--b','linewidth',linewidth)
xlabel('time [s]','fontsize',fontsize); ylabel('u [m]','fontsize', fontsize)
title(['normalized plotted dof location: ',num2str(plotDofLocation(1)),',',...
    num2str(plotDofLocation(2)),',  dof: ',num2str(dofPl)]);
ax = gca; grid on; ax.FontSize = fontsize; legend('HFM','RED','LIN','fontsize',fontsize)
subplot 312
semilogy(dynFull.nlin.time,100*errROM,'-k','linewidth',linewidth)
title('normalized RMS error');
xlabel('time [s]','fontsize', fontsize);
ylabel('error % ', 'fontsize', fontsize);
ax = gca; grid on; ax.FontSize = fontsize;
text('Units', 'Normalized', 'Position', [0.1, 0.1], 'string',...
    description, 'fontsize', 16)
subplot 313
plot(dynRed.nlin.time,errRes,'-b','linewidth',linewidth)
title('full residual at ROM solution');
xlabel('time [s]','fontsize', fontsize);
ylabel('error', 'fontsize', fontsize);
ax = gca; grid on; ax.FontSize = fontsize;

%error divided in membrane dofs, rotational dofs and transversal dofs
figure('units','normalized','position',[0.3,0.1,0.6,0.7]); 
subplot 311
hold on
semilogy(dynFull.nlin.time,100*errROMwDisp,'-k','linewidth',linewidth)
title('normalized RMS out of plane disp. error');
xlabel('time [s]','fontsize', fontsize);
ylabel('error % ', 'fontsize', fontsize);
ax = gca; grid on; ax.FontSize = fontsize;
text('Units', 'Normalized', 'Position', [0.1, 0.1], 'string',...
    description, 'fontsize', 16)
subplot 312
hold on
semilogy(dynFull.nlin.time,100*errROMmembrane,'-k','linewidth',linewidth)
title('normalized RMS membrane error');
xlabel('time [s]','fontsize', fontsize);
ylabel('error % ', 'fontsize', fontsize);
ax = gca; grid on; ax.FontSize = fontsize;
subplot 313
hold on
semilogy(dynFull.nlin.time,100*errROMrotation,'-k','linewidth',linewidth)
title('normalized RMS rotation error');
xlabel('time [s]','fontsize', fontsize);
ylabel('error % ', 'fontsize', fontsize);
ax = gca; grid on; ax.FontSize = fontsize;


%% bending moment history
%extract internal bending action at given node
stressNodeLocation = [0.3,0.2];
nodeSetMoment = find_node(stressNodeLocation(1)*Lx,stressNodeLocation(2)*Ly,[],nodes);

method = 'bending_moments';
nOut = 3;

%extract bending moment field from Full Solution
timeSamplsInt = dynFull.nlin.time;
nTimeSamplesInt = length(timeSamplsInt);
TSamplsInt = zeros(myMesh.nNodes,nTimeSamplesInt);
gradTSamplsInt = zeros(myMesh.nNodes,nTimeSamplesInt);
for ii = 1:nTimeSamplesInt
    TSamplsInt(:,ii) = Tdynt(timeSamplsInt(ii));
    gradTSamplsInt(:,ii) = gradTt(timeSamplsInt(ii));
end

bendingMomentFull = Assembly.get_field(elements,nodeSetMoment,method,nOut,uDynFull,TSamplsInt,gradTSamplsInt);

%extract bending moment field from Red Solution
timeSamplsInt = dynRed.nlin.time;
nTimeSamplesInt = length(timeSamplsInt);
TSamplsInt = zeros(myMesh.nNodes,nTimeSamplesInt);
gradTSamplsInt = zeros(myMesh.nNodes,nTimeSamplesInt);
for ii = 1:nTimeSamplesInt
    TSamplsInt(:,ii) = Tdynt(timeSamplsInt(ii));
    gradTSamplsInt(:,ii) = gradTt(timeSamplsInt(ii));
end

bendingMomentRed = Assembly.get_field(elements,nodeSetMoment,method,nOut,uDynRed,TSamplsInt,gradTSamplsInt);

%% plot comparison (bending moment)

figure('units','normalized','position',[0.3,0.1,0.6,0.7]); hold on
plot(dynFull.nlin.time,squeeze(bendingMomentFull(1,1,:)),'-k','linewidth',linewidth)
plot(dynRed.nlin.time,squeeze(bendingMomentRed(1,1,:)),'--r','linewidth',linewidth)
xlabel('time [s]','fontsize',fontsize); ylabel('M/l [N]','fontsize', fontsize)
title(['M_x, normalized plotted dof location: ',num2str(stressNodeLocation(1)),',',...
    num2str(stressNodeLocation(2))]);
ax = gca; grid on; ax.FontSize = fontsize; legend('HFM','RED','fontsize',fontsize);

%% Compute error on bending moments
method = 'bending_moments';
nodeSetMoment = 'all'; %compute error on bending moment over all the nodes
nOut = 3;

%extract bending moment field from Full Solution
timeSamplsInt = dynFull.nlin.time;
nTimeSamplesInt = length(timeSamplsInt);
TSamplsInt = zeros(myMesh.nNodes,nTimeSamplesInt);
gradTSamplsInt = zeros(myMesh.nNodes,nTimeSamplesInt);
for ii = 1:nTimeSamplesInt
    TSamplsInt(:,ii) = Tdynt(timeSamplsInt(ii));
    gradTSamplsInt(:,ii) = gradTt(timeSamplsInt(ii));
end

bendingMomentFull = Assembly.get_field(elements,nodeSetMoment,method,nOut,uDynFull,TSamplsInt,gradTSamplsInt);

%extract bending moment field from Red Solution
timeSamplsInt = dynRed.nlin.time;
nTimeSamplesInt = length(timeSamplsInt);
TSamplsInt = zeros(myMesh.nNodes,nTimeSamplesInt);
gradTSamplsInt = zeros(myMesh.nNodes,nTimeSamplesInt);
for ii = 1:nTimeSamplesInt
    TSamplsInt(:,ii) = Tdynt(timeSamplsInt(ii));
    gradTSamplsInt(:,ii) = gradTt(timeSamplsInt(ii));
end

bendingMomentRed = Assembly.get_field(elements,nodeSetMoment,method,nOut,uDynRed,TSamplsInt,gradTSamplsInt);

errMomentX = RMS_error(squeeze(bendingMomentFull(:,1,:)),dynFull.nlin.time,...
    squeeze(bendingMomentRed(:,1,:)),dynRed.nlin.time);
errMomentY = RMS_error(squeeze(bendingMomentFull(:,2,:)),dynFull.nlin.time,...
    squeeze(bendingMomentRed(:,2,:)),dynRed.nlin.time);

%plot error
fontsize = 15;
linewidth = 2;

figure('units','normalized','position',[0.3,0.1,0.6,0.7]); hold on
subplot 311
semilogy(dynFull.nlin.time,100*errMomentX ,'-k','linewidth',linewidth)
title('normalized RMS error on Moment Mx');
xlabel('time [s]','fontsize', fontsize);
ylabel('error % ', 'fontsize', fontsize);
ax = gca; grid on; ax.FontSize = fontsize;

subplot 312
semilogy(dynFull.nlin.time,100*errMomentY ,'-k','linewidth',linewidth)
title('normalized RMS error on Moment My');
xlabel('time [s]','fontsize', fontsize);
ylabel('error % ', 'fontsize', fontsize);
ax = gca; grid on; ax.FontSize = fontsize;

subplot 313
plot(dynFull.nlin.time,squeeze(bendingMomentFull(:,1,:)),'-k');
hold on
plot(dynRed.nlin.time,squeeze(bendingMomentRed(:,1,:)),'--r');


