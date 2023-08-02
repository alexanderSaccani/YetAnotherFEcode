% author: ALEXANDER SACCANI - PHD CANDIDATE, ETH ZURICH, 
% created: 1.12.2022
% last modified: 30.01.2023

% Shallow arc beam (the one implemented in Shobhit's paper: Jain, Tiso, Model order
% reduction for temperature dependent nonlinear mechanical systems: A
% multiple scales approach, Journal of Sound and Vibration, 2020


clearvars ; close all; clc
addpath('C:\Users\Pc\OneDrive - ETH Zurich\PhD\YetAnotherFEcode\ROM thermoelastic\general purpose functions 1')
%% Thermoelastic clamped plate

%Material__________________________________________________________________
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
Assembly = Assembly(myMesh,true);
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

%%
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

% VMs = cell(nVMs,1);
% for ii = 1:nVMs  
%    VMs{ii,1} = reshape(VMsC(:,ii),5,[]).';    
% end
% 
% for ii = 1:nVMs
%     figure; hold on;
%     VM2plot = VMs{ii};
%     PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [VM2plot(:,1),...
%         VM2plot(:,2), VM2plot(:,3)], 'factor', 15*max(nodes(:,2)));
%     title(['VM ',num2str(ii), ' f: ', num2str(fcold(ii)),' Hz']);
% end
% 
% nMDs = size(MDsC,2);
% MDs = cell(nMDs,1);
% for ii = 1:nMDs  
%    MDs{ii,1} = reshape(MDsC(:,ii),5,[]).';    
% end
% 
% for ii = 1:nMDs
%     figure; hold on;
%     MD2plot = MDs{ii};
%     PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [MD2plot(:,1),...
%         MD2plot(:,2), MD2plot(:,3)], 'factor',1*max(nodes(:,2)));
%     title(['MD ',num2str(MDsNames(ii,1)),num2str(MDsNames(ii,2))]);
% end


%cold basis of VMs and MDs
Vc = orth(Assembly.constrain_vector([VMsC,MDsC]));

%% Define shape of temperature distribution

%define the shape of T distribution
sigma3X = 0.5*Lx;
sigma3Y = 0.5*Lx;
theta = 20; %inclination w.r.t. x axis

Tmax = 110; %peak value of T

%T as function of p
T_p = @(p) T_p_vararg(p,sigma3X,sigma3Y,theta,nodes,Tmax);

%plot Tshape
xc = Lx/3; %center of pulse
yc = Ly/3; %center of pulse
TshapePlot = T_p([xc,yc]);

figure
hf = PlotFieldonMesh(nodes,elementsPlot,TshapePlot);
title('Temperature distribution shape')


%% Hot modes
pStar = [0.4*Lx,0.4*Ly]';
T_samples = T_p(pStar);
gradT_samples = zeros(myMesh.nNodes,1);

nVMsHot = 4;

logic_MD = 1;
VMMD = VMMD_thermal(Assembly,nVMsHot,logic_MD,T_samples,gradT_samples);
fcold = VMMD.omega/2/pi;
VMsH = VMMD.VMs{1};
MDsH = VMMD.MDs.data{1};
MDsNames = VMMD.MDs.names;

VMsHot = cell(nVMsHot,1);
for ii = 1:nVMsHot  
   VMsHot{ii,1} = reshape(VMsH(:,ii),5,[]).';    
end

VMlist = [1,2,3,4];
for ii = 1:length(VMlist)
    figure; hold on;
    VM2plot = VMsHot{VMlist(ii)};
    PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [VM2plot(:,1),...
        VM2plot(:,2), VM2plot(:,3)], 'factor', 15*max(nodes(:,2)));
    title(['VM ',num2str(VMlist(ii)), ' f: ', num2str(fcold(VMlist(ii))),' Hz']);
end

% nMDs = size(MDsC,2);
% MDs = cell(nMDs,1);
% for ii = 1:nMDs  
%    MDs{ii,1} = reshape(MDsC(:,ii),5,[]).';    
% end
% 
% for ii = 1:nMDs
%     figure; hold on;
%     MD2plot = MDs{ii};
%     PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [MD2plot(:,1),...
%         MD2plot(:,2), MD2plot(:,3)], 'factor',1*max(nodes(:,2)));
%     title(['MD ',num2str(MDsNames(ii,1)),num2str(MDsNames(ii,2))]);
% end

% %% Define parameter grid (xc,yc)
% % %rectangular sampling
% % gridSamplesX = 4;
% % gridSamplesY = 4;
% % 
% % delX = Lx/(gridSamplesX-1);
% % delY = Ly/(gridSamplesY-1);
% % 
% % xcLocations = 0:delX:delX*(gridSamplesX-1);
% % ycLocations = 0:delY:delY*(gridSamplesY-1);
% % 
% % [X,Y] = meshgrid(xcLocations,ycLocations);
% % p = [reshape(X,[],1),reshape(Y,[],1)]';
% 
% %diagonal sampling
% numberTrainingSamples = 20;
% p1 = 0:Lx/(numberTrainingSamples-1):Lx;%xc
% p2 = 0:Ly/(numberTrainingSamples-1):Ly;%yc
% p = [p1;p2];
% p = [[-11*Lx/(numberTrainingSamples-1);-11*Ly/(numberTrainingSamples-1)],...
%     [-10*Lx/(numberTrainingSamples-1);-10*Ly/(numberTrainingSamples-1)],...
%     [-9*Lx/(numberTrainingSamples-1);-9*Ly/(numberTrainingSamples-1)],...
%     [-8*Lx/(numberTrainingSamples-1);-8*Ly/(numberTrainingSamples-1)],...
%     [-7*Lx/(numberTrainingSamples-1);-7*Ly/(numberTrainingSamples-1)],...
%     [-6*Lx/(numberTrainingSamples-1);-6*Ly/(numberTrainingSamples-1)],...
%     [-5*Lx/(numberTrainingSamples-1);-5*Ly/(numberTrainingSamples-1)],...
%     [-4*Lx/(numberTrainingSamples-1);-4*Ly/(numberTrainingSamples-1)],...
%     [-3*Lx/(numberTrainingSamples-1);-3*Ly/(numberTrainingSamples-1)],...
%     [-2*Lx/(numberTrainingSamples-1);-2*Ly/(numberTrainingSamples-1)],...
% [-Lx/(numberTrainingSamples-1);-Ly/(numberTrainingSamples-1)],p];
% 
% figure
% plot(p(1,:),p(2,:),'x','markersize',10,'linewidth',2)
% hold on
% plot([0,Lx,Lx,0,0],[0,0,Ly,Ly,0],'-b','linewidth',2);
% xlabel('[m]'); ylabel('[m]');
% title('training grid')

%%
%diagonal sampling
numberTrainingSamples = 40;
p1 = 0:Lx/(numberTrainingSamples-1):Lx;%xc
p2 = 0:Ly/(numberTrainingSamples-1):Ly;%yc
p = [p1;p2];
p = [[-11*2*Lx/(numberTrainingSamples-1);-11*2*Ly/(numberTrainingSamples-1)],...
    [-10*2*Lx/(numberTrainingSamples-1);-10*2*Ly/(numberTrainingSamples-1)],...
    [-9*2*Lx/(numberTrainingSamples-1);-9*2*Ly/(numberTrainingSamples-1)],...
    [-8*2*Lx/(numberTrainingSamples-1);-8*2*Ly/(numberTrainingSamples-1)],...
    [-7*2*Lx/(numberTrainingSamples-1);-7*2*Ly/(numberTrainingSamples-1)],...
    [-6*2*Lx/(numberTrainingSamples-1);-6*2*Ly/(numberTrainingSamples-1)],...
    [-5*2*Lx/(numberTrainingSamples-1);-5*2*Ly/(numberTrainingSamples-1)],...
    [-4*2*Lx/(numberTrainingSamples-1);-4*2*Ly/(numberTrainingSamples-1)],...
    [-3*2*Lx/(numberTrainingSamples-1);-3*2*Ly/(numberTrainingSamples-1)],...
    [-2*2*Lx/(numberTrainingSamples-1);-2*2*Ly/(numberTrainingSamples-1)],...
[-Lx/(numberTrainingSamples-1);-Ly/(numberTrainingSamples-1)],p];

figure
plot(p(1,:),p(2,:),'x','markersize',10,'linewidth',2)
hold on
plot([0,Lx,Lx,0,0],[0,0,Ly,Ly,0],'-b','linewidth',2);
xlabel('[m]'); ylabel('[m]');
title('training grid')

%% Construct Multiple ROMs
nVMs = 6;
thermEq = true;
orthog = true;
reordBasis = true;
gradFun.flag = true; %must put it true because it is shell element
gradFun.fun = @(dummy) zeros(Assembly.Mesh.nNodes,1);
gradFun.param = zeros(1,length(p));
unconstrainBasis = true;
refBasis = 1;
[ROMs.models,sampledVMs,sampledMDs,sampledEqDisp,natFreq] = multiple_ROMs_pthermal(Assembly,...
            p, T_p, gradFun, nVMs, thermEq, orthog, reordBasis, unconstrainBasis, refBasis);
  
%do this inside of function next time
ROMs.psamples = p;

%% Plot the thermal equilibrium
eq = sampledEqDisp{23};
eq = reshape(eq,5,[]).';
figure
PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [eq(:,1),...
            eq(:,2), eq(:,3)], 'factor',10);

%%
%% Test on thermal trajectory

%FORCING___________________________________________________________________
%shape
%pressure (triangular)
nodesX = nodes(:,1);
nodesY = nodes(:,2);
FmaxX = 900/length(nodes);
FmaxY = 900/length(nodes);
Fmin = 0;
Fz = FmaxX-nodesX*FmaxX/Lx + FmaxY - nodesY*FmaxY/Ly;
F0 = zeros(myMesh.nDOFs,1);
F0(3:5:end) = Fz;

%plot forcing shape
figure
hf = PlotFieldonMesh(nodes,elementsPlot,Fz);
title('Forcing distribution shape (out of plane)')

fForc = 150; %close to second frequency
omForc = fForc*2*pi;
TForc = 1/fForc;

Fext = @(t) 0*F0*cos(omForc*t);
%%
%TEMPERATURE VARIATION_____________________________________________________
thetaTraj = 180/pi*atan(Ly/Lx);%30; %inclination of center trajectory
h0 = 0;%Ly/4;
time2Tranverse = 200*TForc; %slow temperature distribution
hoffset = 0.15;%Ly/3;

p_fun_t = @(t)p_fun_t_vararg(thetaTraj,h0,time2Tranverse,Ly,hoffset,t);

figure
plot(p(1,:),p(2,:),'x','markersize',10,'linewidth',2)
hold on
plot([0,Lx,Lx,0,0],[0,0,Ly,Ly,0],'-b','linewidth',2);
xlabel('[m]'); ylabel('[m]');
title('training grid')
timePlot = linspace(0,time2Tranverse,1000);
pHist = p_fun_t(timePlot); hold on;
plot(pHist(1,:),pHist(2,:),'-r','linewidth',2)
%PlotMesh([nodes,zNodes],elementsPlot);
axis([-hoffset/tan(thetaTraj*pi/180),Lx,-hoffset,Ly]);

Tdynt = @(t)T_p(p_fun_t(t));
gradTt = @(t) zeros(Assembly.Mesh.nDOFs,1);

%INTEGRATION TIME__________________________________________________________
tint = 0.7*time2Tranverse;


%plot dynamic T over the analysis
tPlot = linspace(0,tint,100);
figure
for ii = 1:length(tPlot)
    tii = tPlot(ii);
    hf = PlotFieldonMesh(nodes,elementsPlot,Tdynt(tii));
    title('Temperature distribution in time')
    clim([0, Tmax]); 
    txt = text('Units', 'Normalized', 'Position', [0.5, 0.9], 'string',...
    ['t: ',num2str(round(tii,5)),'s'], 'fontsize', 12); drawnow limitrate
    delete(txt)
end
txt = text('Units', 'Normalized', 'Position', [0.4, 0.9], 'string',...
    ['t: ',num2str(round(tii,5)),'s'], 'fontsize', 12);


%% RUN HFM
%HFM_______________________________________________________________________

runHFM = false;

if runHFM 

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
tic
TI_NL.Integrate(q0,qd0,qdd0,tint,residualFull);
toc

% project to full space
uDynFull = TI_NL.Solution.q;

%save solution
uDynFull = Assembly.unconstrain_vector(uDynFull); 
dynFull.nlin.disp = decodeDofsNodes(uDynFull,nNodes,5); % (node, dof of node, tsamp)
dynFull.nlin.time = TI_NL.Solution.time;

save('saved_simulations/HFRun_noForc.mat','dynFull');

else
    
filenameHFRun = 'saved_simulations/HFRun.mat';
load(filenameHFRun)

end

%% Simulation with interpolated basis

RUNROM = true;
if RUNROM

    % Initial condition: equilibrium
    nDofRed = size(ROMs.models{1}.V,2);

    q0 = zeros(nDofRed,1);
    qd0 = zeros(nDofRed,1);
    qdd0 = zeros(nDofRed,1);

    % Instantiate object for nonlinear time integration
    TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);


    % % Constant basis
    % for ii = 1:numberTrainingSamples
    %
    % ROMs.models{ii}.V = ROMs.models{2}.V;
    %
    % end

    %prepare array for interpolation
    nROMs = length(ROMs.models);
    nVecBasis = size(ROMs.models{1}.V,2);
    nDofsFull = size(ROMs.models{1}.V,1);
    array4Int = zeros(nROMs,nDofsFull,nVecBasis);
    for ii = 1:nROMs
        array4Int(ii,:,:) = ROMs.models{ii}.V;
    end

    % nonlinear Residual evaluation function handle
    residualInterpBasis = @(q,qd,qdd,t)thermal_residual_interpolation_mod(q,qd,qdd,t,Assembly,ROMs,p_fun_t,Fext,array4Int,Tdynt,gradTt);

    % integrate equations with Newmark scheme
    tic
    TI_NL.Integrate(q0,qd0,qdd0,tint,residualInterpBasis);
    toc


    % project to full space
    uDynRedInt = TI_NL.Solution.q;
    timeRedInt = TI_NL.Solution.time;

    uDynRedInt = to_full_int_mod(uDynRedInt ,timeRedInt,p_fun_t,ROMs,array4Int);

    %save solution
    dynRedInt.nlin.disp = decodeDofsNodes(uDynRedInt,nNodes,5); % (node, dof of node, tsamp)
    dynRedInt.nlin.time = timeRedInt;


    save('saved_simulations/RedRun_noForc_1vec.mat','dynRedInt');


else

    filenameRedRun = 'saved_simulations/RedRun.mat';
    load(filenameRedRun)

end

%% quasistatic solution from interpolation of equilibria

uQuasist = interpolated_quasistatic(timeRedInt,ROMs.psamples,sampledEqDisp,p_fun_t);

dynQuasist.nlin.disp = decodeDofsNodes(uQuasist,nNodes,5); % (node, dof of node, tsamp)
dynQuasist.nlin.time = timeRedInt;


%% Plot results
%displacements
plotDofLocation = [0.2,0.2]; %percentage of plate length (along x and y)
node2plot = find_node(plotDofLocation(1)*Lx,plotDofLocation(2)*Ly,[],nodes); % node to plot
dofPl = 1;

fontsize = 15;
linewidth = 2;

figure('units','normalized','position',[0.3,0.1,0.6,0.7]); hold on
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,dofPl,:)),'-k','linewidth',linewidth)
plot(dynRedInt.nlin.time,squeeze(dynRedInt.nlin.disp(node2plot,dofPl,:)),'--y','linewidth',linewidth)
plot(dynRedInt.nlin.time,squeeze(dynQuasist.nlin.disp(node2plot,dofPl,:)),'--m','linewidth',linewidth)
%plot(dynRedGlobal.nlin.time,squeeze(dynRedGlobal.nlin.disp(node2plot,dofPl,:)),'-g','linewidth',1.5)
xlabel('time [s]','fontsize',fontsize); ylabel('u [m]','fontsize', fontsize)
title(['normalized plotted dof location: ',num2str(plotDofLocation(1)),',',...
    num2str(plotDofLocation(2)),',  dof: ',num2str(dofPl)]);
ax = gca; grid on; ax.FontSize = fontsize; legend('HFM','RED int','quasistInt');%,'RED global','fontsize',fontsize);




%% Temperature function
%single gaussian
function T = gaussian_temperature_profile(sigma3X,sigma3Y,thetaDeg,xc,yc,nodes,Tmax)

%dependent parameteres
sigmaX = sigma3X/3;
sigmaY = sigma3Y/3;
theta = pi*thetaDeg/180;

a = cos(theta)^2/(2*sigmaX^2) + sin(theta)^2/(2*sigmaY^2);
b = - sin(2*theta)/(4*sigmaX^2) + sin(2*theta)/(4*sigmaY^2);
c = sin(theta)^2/(2*sigmaX^2) + cos(theta)^2/(2*sigmaY^2);

nodesX = nodes(:,1);
nodesY = nodes(:,2);

T =  Tmax*exp( -( a*(nodesX-xc).^2 +2*b*(nodesX-xc).*(nodesY-yc) +c*(nodesY-yc).^2) ); %gaussian with center at xc,yc inclined of angle theta w.r.t. x axis

end

%Tshape of p (parameter)
function T = T_p_vararg(p,sigma3X,sigma3Y,theta,nodes,Tmax)

xc = p(1);
yc = p(2);

T = gaussian_temperature_profile(sigma3X,sigma3Y,theta,xc,yc,nodes,Tmax);

end

%temporal variation of p___________________________________________________

%variation of xc, yc on a line
function [xc,yc] =  linear_center_variation(theta,h0,time2Tranverse,Ly,hoffset,t)

yc = (Ly-h0+hoffset)*t/time2Tranverse + h0-hoffset;
xc = 1/(tan(theta*pi/180))*(yc-h0);

end

function p = p_fun_t_vararg(theta,h0,time2Tranverse,Ly,hoffset,t)

[xc,yc] = linear_center_variation(theta,h0,time2Tranverse,Ly,hoffset,t);
p = [xc;yc];

end





