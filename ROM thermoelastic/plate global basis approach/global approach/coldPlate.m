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
height_midspan = 8e-3; %height of midspan 

%Mesh______________________________________________________________________
nx = 15; % # elements along x
ny = 8; % # elements along y

eltype = 'QUAD8'; %'QUAD4'
myElementConstructor = @()ThermQuad8Shell(thickness, myMaterial,2);
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
%myMesh.set_essential_boundary_condition([nset{2}],1:3,0)
myMesh.set_essential_boundary_condition([nset{3}],1:3,0)
%myMesh.set_essential_boundary_condition([nset{4}],1:3,0)

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


%% VMs and MDs Cold

T_samples = zeros(myMesh.nNodes,1);
gradT_samples = zeros(myMesh.nNodes,1);

nVMs = 4;

logic_MD = 1;
VMMD = VMMD_thermal(Assembly,nVMs,logic_MD,T_samples,gradT_samples);
fcold = VMMD.omega/2/pi;
VMsC = VMMD.VMs{1};
MDsC = VMMD.MDs.data{1};
MDsNames = VMMD.MDs.names;
eqPoint = VMMD.eq;

VMs = cell(nVMs,1);
for ii = 1:nVMs  
   VMs{ii,1} = reshape(VMsC(:,ii),5,[]).';    
end

for ii = 1:nVMs
    figure; hold on;
    VM2plot = VMs{ii};
    PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [VM2plot(:,1),...
        VM2plot(:,2), VM2plot(:,3)], 'factor', 15*max(nodes(:,2)));
    title(['VM ',num2str(ii), ' f: ', num2str(fcold(ii)),' Hz']);
end

nMDs = size(MDsC,2);
MDs = cell(nMDs,1);
for ii = 1:nMDs  
   MDs{ii,1} = reshape(MDsC(:,ii),5,[]).';    
end

for ii = 1:nMDs
    figure; hold on;
    MD2plot = MDs{ii};
    PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [MD2plot(:,1),...
        MD2plot(:,2), MD2plot(:,3)], 'factor',1*max(nodes(:,2)));
    title(['MD ',num2str(MDsNames(ii,1)),num2str(MDsNames(ii,2))]);
end


%cold basis of VMs and MDs
Vc = orth(Assembly.constrain_vector([VMsC,MDsC]));


disp(['Cold basis has size ',num2str(size(Vc,2))])

%% Test on thermal trajectory

%FORCING___________________________________________________________________
%shape
%pressure (triangular)
nodesX = nodes(:,1);
nodesY = nodes(:,2);
FmaxX = 3000/length(nodes);
FmaxY = 3000/length(nodes);
Fmin = 0;
Fz = FmaxX-nodesX*FmaxX/Lx + FmaxY - nodesY*FmaxY/Ly;
F0 = zeros(myMesh.nDOFs,1);
F0(3:5:end) = Fz;

%plot forcing shape
figure
hf = PlotFieldonMesh(nodes,elementsPlot,Fz);
title('Forcing distribution shape (out of plane)')

fForc = 140; %close to second frequency
omForc = fForc*2*pi;
TForc = 1/fForc;

Fext = @(t) F0*cos(omForc*t);

Tdynt = @(t)zeros(Assembly.Mesh.nNodes,1);
gradTt = @(t) zeros(Assembly.Mesh.nNodes,1);

%INTEGRATION TIME__________________________________________________________
tint = 10*TForc;
%HFM_______________________________________________________________________

runHFM = true;

if runHFM 

% settings for integration
h = 2*pi/omForc/35; % time step for integration

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
tic
TI_NL.Integrate(q0,qd0,qdd0,tint,residualFull);
toc

% project to full space
uDynFull = TI_NL.Solution.q;

%save solution
uDynFull = Assembly.unconstrain_vector(uDynFull); 
dynFull.nlin.disp = decodeDofsNodes(uDynFull,nNodes,5); % (node, dof of node, tsamp)
dynFull.nlin.time = TI_NL.Solution.time;

save('saved_simulations/HFRunColdPlate.mat','dynFull');

else
    
filenameHFRun = 'HFRunColdPlate.mat';
load(filenameHFRun)

end


%%
%plot HFM simulation
dof2plot = 3;
plot_dof_location = [0.3,0.3]; %percentage of length of beam
node2plot = find_node(plot_dof_location(1)*Lx,plot_dof_location(2)*Ly,[],nodes); % node to plot
figure
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,dof2plot,:)))
hold on
%%
%LINEARIZED HFM____________________________________________________________
% %% integrate linearized model
% % settings for integration
% h = 2*pi/omForc/25; % time step for integration
% 
% % Initial condition: equilibrium
% nDofFull = length(Assembly.Mesh.EBC.unconstrainedDOFs);
% 
% q0 = zeros(nDofFull,1);
% qd0 = zeros(nDofFull,1);
% qdd0 = zeros(nDofFull,1);
% 
% % Instantiate object for nonlinear time integration
% TI_LIN = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
% 
% % nonlinear Residual evaluation function handle
% residualFullLin = @(q,qd,qdd,t)residual_thermal_linear(q,qd,qdd,t,Assembly,Fext,Tdynt,gradTt);
% 
% % integrate equations with Newmark scheme
% TI_LIN.Integrate(q0,qd0,qdd0,tint,residualFullLin);
% 
% % project to full space
% uDynFullLin = TI_LIN.Solution.q;
% 
% %save solution
% uDynFullLin = Assembly.unconstrain_vector(uDynFullLin); 
% dynFull.lin.disp = decodeDofsNodes(uDynFullLin,nNodes,5); % (node, dof of node, tsamp)
% dynFull.lin.time = TI_LIN.Solution.time;
% 
% plot_dof_location = [0.3,0.3]; %percentage of length of beam
% node2plot = find_node(plot_dof_location(1)*Lx,plot_dof_location(2)*Ly,[],nodes); % node to plot
% 
% figure
% plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,3,:)))
%  hold on
%  plot(dynFull.lin.time,squeeze(dynFull.lin.disp(node2plot,3,:)))
%  axis([0,tint,-0.02,0.02])
% % legend('NL','LIN');


%% Integrate nonlinear ROM (cold basis)
% ROM (Cold Basis)
V = Vc;

V = Assembly.unconstrain_vector(V);
redAssembly = ReducedAssembly(myMesh,V); %reduced assembly
redAssembly.DATA.M = V.'*Assembly.DATA.M*V;
redAssembly.DATA.C = V.'*Assembly.DATA.C*V;

%INTEGRATION________
runROM = true;

if runROM 

% settings for integration
h = 2*pi/omForc/35; % time step for integration

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
dynRedCold.nlin.disp = decodeDofsNodes(uDynRed,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dynRedCold.nlin.time = TI_NL_r.Solution.time;

save('saved_simulations/ROMColdColdPlate.mat','dynRedCold');

else
filenameROMRun = 'ROMColdColdPlate.mat';
load(filenameROMRun)
end

%%
%displacements
plotDofLocation = [0.4,0.4]; %percentage of plate length (along x and y)
node2plot = find_node(plotDofLocation(1)*Lx,plotDofLocation(2)*Ly,[],nodes); % node to plot
dofPl = 3;

%description = 'reduction basis: first 5 cold VMs and MDs + static T + static TF, nDofRed = 22 ';
fontsize = 15;
linewidth = 2;

figure('units','normalized','position',[0.3,0.1,0.6,0.7]); hold on
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,dofPl,:)),'-k','linewidth',linewidth)
plot(dynRedCold.nlin.time,squeeze(dynRedCold.nlin.disp(node2plot,dofPl,:)),'--r','linewidth',linewidth)
xlabel('time [s]','fontsize',fontsize); ylabel('u [m]','fontsize', fontsize)
title(['normalized plotted dof location: ',num2str(plotDofLocation(1)),',',...
    num2str(plotDofLocation(2)),',  dof: ',num2str(dofPl)]);
ax = gca; grid on; ax.FontSize = fontsize; legend('HFM','RED','fontsize',fontsize);

