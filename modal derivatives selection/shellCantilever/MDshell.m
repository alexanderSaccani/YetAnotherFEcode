close all
clearvars
clc
%% Construct model

%Material _________________________________________________________________
E       = 71e9;     % Young's modulus [Pa]
rho     = 2795;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 

myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);

%Geometry__________________________________________________________________
% Lx = 0.65;
% Ly = .125;
% thickness = .003;     % thickness of plate
% height_midspan = 5e-10; %height of midspan 

Lx = 0.55;
Ly = .06;
thickness = .00225;     % thickness of plate
height_midspan = 5e-10; %height of midspan 

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
elementsPlot = elements(:,1:4); %for plots
PlotMesh([nodes,zNodes],elementsPlot);

%Boundary conditions: all sides are clamped________________________________
myMesh.set_essential_boundary_condition([nset{1}],1:5,0)

%Assembly__________________________________________________________________
Assembly = Assembly(myMesh);
M = Assembly.mass_matrix();
u0 = zeros( myMesh.nDOFs, 1);

%Compute and store Mass matrix, K matrix
[K,~] = Assembly.tangent_stiffness_and_force(u0);

% store matrices
Assembly.DATA.K = K; %K
Assembly.DATA.M = M; %mass matrix


%% VMs analysis 

nVMs = 5; % first n_VMs modes with lowest frequency calculated 

Kconstr = Assembly.constrain_matrix(Assembly.DATA.K);
Mconstr = Assembly.constrain_matrix(Assembly.DATA.M);

[VMs,om] = eigs(Kconstr, Mconstr, nVMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);

VMs = VMs(:,ind);
for ii = 1:nVMs
    VMs(:,ii) = VMs(:,ii)/max(sqrt(sum(VMs(:,ii).^2,2)));
end
VMs = Assembly.unconstrain_vector(VMs);

%plot VMs 
VMsPlot = cell(nVMs,1); 
for ii = 1:nVMs
    VMsPlot{ii,1} = reshape(VMs(:,ii),5,[]).';
end

VMs2Plot = 1:5; 
for ii = 1:length(VMs2Plot)
    figure; hold on; VM2plotii = VMsPlot{VMs2Plot(ii)};
    PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot,[VM2plotii(:,1),...
        VM2plotii(:,2), VM2plotii(:,3)], 'factor',1.5*max(nodes(:,2)));
    title(['VM ',num2str(VMs2Plot(ii)),', f: ',num2str(f0(VMs2Plot(ii))),' Hz']);
end

%% Rayleigh Damping
r1 = 0.01; %damping of first mode
r2 = 0.01;  %damping of second mode

r = [r1;r2];
om1 = 2*pi*f0(1);
om2 = 2*pi*f0(2);

dampParam = [1 om1^2; 1 om2^2]\r;

C = dampParam(1)*M + dampParam(2)*K;

Assembly.DATA.C = C;


%% Static Analysis
% Nodal force
F = zeros(myMesh.nDOFs,1);
nf = find_node(0.4*Lx,0.25*Ly,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F(node_force_dofs(3)) = 75;


[ulin,unl] = static_equilibrium(Assembly, F, 'display', 'iter-detailed');
UNL = reshape(unl,5,[]).';        % Nonlinear response
ULIN = reshape(ulin,5,[]).';	% Linear response


% PLOT
elementPlot = elements(:,1:4);
figure('units','normalized','position',[.2 .1 .6 .8])
scale = 1;
PlotMesh(nodes, elementPlot, 0);
PlotFieldonDeformedMesh([nodes,zNodes],elementPlot,[ULIN(:,1),ULIN(:,2),ULIN(:,3)],'factor',scale,'color','k');
hold on
PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot,[UNL(:,1),UNL(:,2),UNL(:,3)],'factor',scale,'color','k');
colormap jet
title(['LINEAR and NONLINEAR STATIC RESPONSE, scale factor: ' num2str(scale) 'x)']);

%% Modal Derivatives

%Static Modal Derivatives
%how many VMs
nVMs4MDs = 4;
[SMDs, names] = modal_derivatives(Assembly, elements, VMs(:,1:nVMs4MDs));

%% Construct ROM
%construct ROM basis
V = orth([VMs,SMDs]); 
nDOFred = size(V,2);

%compute tensors
tensors = reduced_tensors_ROM(Assembly, elements, V);

%create object reduced assembly
RedAssembly = ReducedAssembly(myMesh,V);
RedAssembly.DATA.M = V'*Assembly.DATA.M*V;
RedAssembly.DATA.C = V'*Assembly.DATA.C*V;

%% Dynamic Analysis
%nodal force
F1 = zeros(myMesh.nDOFs,1);
nf = find_node(0.4*Lx,0.25*Ly,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F1(node_force_dofs(3)) = 10;

%frequency content
fForc1 = 34; %excite second mode 
omForc1 = 2*pi*fForc1;
TForc1 = 1/fForc1;

Fext = @(t) F1*cos(omForc1*t);

%integration time
tint = 6*TForc1;

hnlin = TForc1/100;

%% Integrate full model nonlinear

% Initial condition: equilibrium
nDofFull = length(Assembly.Mesh.EBC.unconstrainedDOFs);

q0 = zeros(nDofFull,1);
qd0 = zeros(nDofFull,1);
qdd0 = zeros(nDofFull,1);

% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',hnlin,'alpha',0.005);

% nonlinear Residual evaluation function handle
residualFull = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,Assembly,Fext);

% integrate equations with Newmark scheme
TI_NL.Integrate(q0,qd0,qdd0,tint,residualFull);

% project to full space
uDynFull = TI_NL.Solution.q;

%save solution
uDynFull = Assembly.unconstrain_vector(uDynFull); 
dynFull.nlin.disp = decodeDofsNodes(uDynFull,nNodes,5); 
dynFull.nlin.time = TI_NL.Solution.time;

%plot solution
plot_dof_location = [1,1]; %percentage of length of beam
node2plot = find_node(plot_dof_location(1)*Lx,plot_dof_location(2)*Ly,[],nodes); % node to plot
figure
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,3,:)))


% tTotAnimation = 10;
% elementPlot = elements(:,1:4);
% frameRate = length(dynFull.nlin.time)/tTotAnimation;
% 
% % figure
% AnimateFieldonDeformedMesh([nodes,zNodes],elementsPlot,...
%  [uDynFull(:,1),uDynFull(:,2),uDynFull(:,3)],'framerate',round(frameRate),'factor',1)

%% integrate linearized model
% settings for integration
hlin = 2*hnlin;

% Initial condition: equilibrium
nDofFull = length(Assembly.Mesh.EBC.unconstrainedDOFs);

q0 = zeros(nDofFull,1);
qd0 = zeros(nDofFull,1);
qdd0 = zeros(nDofFull,1);

% Instantiate object for nonlinear time integration
TI_LIN = ImplicitNewmark('timestep',hlin,'alpha',0.005,'linear',true);

% nonlinear Residual evaluation function handle
residualFullLin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,Assembly,Fext);

% integrate equations with Newmark scheme
TI_LIN.Integrate(q0,qd0,qdd0,tint,residualFullLin);

% project to full space
uDynFullLin = TI_LIN.Solution.q;

%save solution
uDynFullLin = Assembly.unconstrain_vector(uDynFullLin); 
dynFull.lin.disp = decodeDofsNodes(uDynFullLin,nNodes,5); % (node, dof of node, tsamp)
dynFull.lin.time = TI_LIN.Solution.time;

plot_dof_location = [1,1]; %percentage of length of beam
node2plot = find_node(plot_dof_location(1)*Lx,plot_dof_location(2)*Ly,[],nodes); % node to plot

figure
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,3,:)))
hold on
plot(dynFull.lin.time,squeeze(dynFull.lin.disp(node2plot,3,:)))
legend('NL','LIN')

%% ROM

Vtest{1} = V;

% ROM
F1red = V'*F1;
FextRed = @(t) F1red*cos(omForc1*t); %red force


dynRed = cell(length(Vtest),1);
errROM{ii} = cell(length(Vtest),1);

for ii = 1:length(Vtest)
    
    V = Vtest{ii}; 
    
    nDofRed = size(V,2);

    q0 = zeros(nDofRed,1);
    qd0 = zeros(nDofRed,1);
    qdd0 = zeros(nDofRed,1);
    
    % Instantiate object for nonlinear time integration
    TI_NL_r = ImplicitNewmark('timestep',hnlin,'alpha',0.005);
    
    % nonlinear Residual evaluation function handle
    residual_r = @(q,qd,qdd,t)residual_nonlinear_ROM_tensor(q,qd,qdd,t,RedAssembly,tensors,FextRed);
    
    % integrate equations with Newmark scheme
    TI_NL_r.Integrate(q0,qd0,qdd0,tint,residual_r);
    
    % project to full space
    uDynRed = V*TI_NL_r.Solution.q;
    vDynRed = V*TI_NL_r.Solution.qd;
    aDynRed = V*TI_NL_r.Solution.qdd;
    
    %save solution
    dynRed{ii}.nlin.disp = decodeDofsNodes(uDynRed,nNodes,5); % (node, dof of node, tsamp)
    dynRed{ii}.nlin.time = TI_NL_r.Solution.time;

    %compute error
    errROM{ii} = RMS_error_ROM(dynFull.nlin.disp,dynFull.nlin.time,dynRed{ii}.nlin.disp,dynRed{ii}.nlin.time,1:2,1:myMesh.nNodes);

end


%% Plot comparison
%displacements
plotDofLocation = [1,1]; %percentage of plate length (along x and y)
node2plot = find_node(plotDofLocation(1)*Lx/2,plotDofLocation(2)*Ly/4,[],nodes); % node to plot
dofPl = 1;

description = '';
fontsize = 15;
linewidth = 2;

figure('units','normalized','position',[0.3,0.1,0.6,0.7]); hold on
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,dofPl,:)),'-k','linewidth',linewidth)
plot(dynRed{1}.nlin.time,squeeze(dynRed{1}.nlin.disp(node2plot,dofPl,:)),'--r','linewidth',linewidth)
plot(dynFull.lin.time,squeeze(dynFull.lin.disp(node2plot,dofPl,:)),'-b')
legend('HFM','ROM','LIN')
