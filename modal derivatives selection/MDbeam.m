% EXAMPLE: beam meshed with 3D element
clear; 
% close all; 
clc

%% construct model

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = .01;     % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain

% Element
myElementConstructor = @()Quad8Element(thickness, myMaterial);

% MESH_____________________________________________________________________
Lx = 0.4;
Ly = .02;
nx = 50;
ny = 2;

[nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,'QUAD8');

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);


myMesh.set_essential_boundary_condition([nset{1}],1:2,0)

% ASSEMBLY ________________________________________________________________
Assembly = Assembly(myMesh);
M = Assembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( myMesh.nDOFs, 1);
[K,~] = Assembly.tangent_stiffness_and_force(u0);

% store matrices
Assembly.DATA.K = K;
Assembly.DATA.M = M;
Assembly.DATA.F0 = zeros(myMesh.nDOFs,1);

                                                  
%% Eigenvalue problem
n_VMs = 30; % first n_VMs modes with lowest frequency calculated 
Kc = Assembly.constrain_matrix(K);
Mc = Assembly.constrain_matrix(M);
[VMs,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
VMs = VMs(:,ind);
for ii = 1:n_VMs
    VMs(:,ii) = VMs(:,ii)/max(sqrt(sum(VMs(:,ii).^2,2)));
end
VMs = Assembly.unconstrain_vector(VMs);

% PLOT
mod2plot = 1:4;
for ii = 1:length(mod2plot)
    mod = mod2plot(ii);
    elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
    figure('units','normalized','position',[.2 .1 .6 .8])
    PlotMesh(nodes, elementPlot, 0);
    v1 = reshape(VMs(:,mod), 2, []).';
    PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
    title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])
end

%% rayleigh damping

r1 = 0.01;
r2 = 0.01;

r = [r1;r2];
om1 = 2*pi*f0(1);
om2 = 2*pi*f0(2);

dampParam = [1 om1^2; 1 om2^2]\r;

C = dampParam(1)*M + dampParam(2)*K;

Assembly.DATA.C = C;
%% Static Modal Derivatives

%how many VMs
nVMs4MDs = 4;
[SMDs, names] = modal_derivatives(Assembly, elements, VMs(:,1:nVMs4MDs));

% mod2plot = 1:3;
% for ii = 1:length(mod2plot)
%     mod = mod2plot(ii);
%     elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
%     figure('units','normalized','position',[.2 .1 .6 .8])
%     PlotMesh(nodes, elementPlot, 0);
%     v1 = reshape(SMDs(:,mod), 2, []).';
%     PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
%     title(['static MD ' num2str(names(mod,1)),num2str(names(mod,2))])
% end

%% Perturbation Modal Derivatives

[PMDs, names, angExcFreq] = dynamic_modal_derivatives(Assembly, elements, VMs(:,1:nVMs4MDs), 2*pi*f0);
PMDs1 = squeeze(PMDs(:,:,1)); %sum of freq
PMDs2 = squeeze(PMDs(:,:,2)); %diff of freq

MDExcFreq = angExcFreq/(2*pi);



% %plot PDM1
% mod2plot = 1:3;
% for ii = 1:length(mod2plot)
%     mod = mod2plot(ii);
%     elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
%     figure('units','normalized','position',[.2 .1 .6 .8])
%     PlotMesh(nodes, elementPlot, 0);
%     v1 = reshape(PMDs1(:,mod), 2, []).';
%     PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
%     title(['dynamic MD1 (sum of freq): ' num2str(names(mod,1)),num2str(names(mod,2))])
% end
% 
% %plot PDM2
% mod2plot = 1:3;
% for ii = 1:length(mod2plot)
%     mod = mod2plot(ii);
%     elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
%     figure('units','normalized','position',[.2 .1 .6 .8])
%     PlotMesh(nodes, elementPlot, 0);
%     v1 = reshape(PMDs2(:,mod), 2, []).';
%     PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
%     title(['Dynamic MD2 (diff of freq): ' num2str(names(mod,1)),num2str(names(mod,2))])
% end

% diff1 = SMDs-PMDs1;
% diff2 = SMDs-PMDs2;
% diff3 = PMDs1 - PMDs2;


%define node set
nodeSet1 = zeros(nx+1,1);
dofSet1 = zeros(nx+1,2);

for ii = 1:nx+1

nii = find_node(Lx/nx*(ii-1),Ly/ny*0,[],myMesh.nodes); % node where to put the force
nodeSet1(ii) = nii;

dofSet1(ii,:) = (get_index(nii, myMesh.nDOFPerNode ))'; % node where to put the force

end

dofxSet1 = dofSet1(:,1);
dofySet1 = dofSet1(:,2);

xnodesSet1 = nodes(nodeSet1,1);
ynodesSet1 = nodes(nodeSet1,2);

MDs2Plot = 1:size(SMDs,2);
% for ii = 1:length(MDs2Plot)
% 
% figure
% subplot 211
% plot(xnodesSet1,SMDs(dofySet1,ii),'-k')
% hold on
% plot(xnodesSet1,PMDs1(dofySet1,ii),'--b')
% plot(xnodesSet1,PMDs2(dofySet1,ii),'-.r')
% title('tranverse direction')
% 
% subplot 212
% plot(xnodesSet1,SMDs(dofxSet1,ii),'-k')
% hold on
% plot(xnodesSet1,PMDs1(dofxSet1,ii),'--b')
% plot(xnodesSet1,PMDs2(dofxSet1,ii),'-.r')
% title('axial direction')
% legend('SMDs','PMDs1: sum Freq','PMDs2: diff Freq')
% 
% freqi = f0(names(ii,1));
% freqj = f0(names(ii,2));
% sumFreq = freqi + freqj;
% diffFreq = abs(freqi - freqj);
% 
% description = ['sum of freqs: ',num2str(sumFreq)];
% text('Units', 'Normalized', 'Position', [0.1, 0.1], 'string',...
%     description, 'fontsize', 10)
% description = ['diff of freqs: ', num2str(diffFreq)];
% text('Units', 'Normalized', 'Position', [0.1, 0.2], 'string',...
%     description, 'fontsize', 10)
% 
% sgtitle(['comparision of MDs: ' num2str(names(ii,1)),num2str(names(ii,2))])
% end

f0inv = 1./f0;
excitedMode = MDExcFreq(:,1)*f0inv(nVMs4MDs+1:end)'; %1 corresponds to sum of frequencies

%% Static Analysis
% Nodal force
F = zeros(myMesh.nDOFs,1);
nf = find_node(Lx,Ly/2,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F(node_force_dofs(2)) = -2500;


u_lin = Assembly.solve_system(K, F);
ULIN = reshape(u_lin,2,[]).';	% Linear response
[ulin,unl] = static_equilibrium(Assembly, F, 'display', 'iter-detailed');
UNL = reshape(unl,2,[]).';        % Nonlinear response


% PLOT
elementPlot = elements(:,1:4);
figure('units','normalized','position',[.2 .1 .6 .8])
scale = 1;%2*Ly/max(abs(UNL(:)));
PlotMesh(nodes, elementPlot, 0);
PlotFieldonDeformedMesh(nodes,elementPlot,UNL,'factor',scale,'color','k');
colormap jet
title(['NONLINEAR STATIC RESPONSE, scale factor: ' num2str(scale) 'x)']);


%% Dynamic Analysis

%shape of load
F1 = zeros(myMesh.nDOFs,1);
nf = find_node(Lx,Ly,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F1(node_force_dofs(2)) = -3800; 

F2 = zeros(myMesh.nDOFs,1);
nf = find_node(Lx,Ly/2,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F2(node_force_dofs(1)) = 0;

%% Analysis of load applied to Linearized System

Ftot = F1 + F2;
om0 = f0(1:20)*2*pi;

OM = linspace(0,max(om0),10000)';
FRF = zeros(length(OM),length(om0));
for ii = 1:length(om0)
    FRF(:,ii)= 1./(-OM.^2+om0(ii)^2+1i*OM*(dampParam(1)+dampParam(2)*om0(ii)^2)); 
end

figure
for ii = 1:length(om0)
   semilogy(OM,abs(FRF(:,ii)),'linewidth',2); hold on;
   grid on; xlabel('\omega'); ylabel('|\eta|');
end
title('modal FRFs of linearized system')

hold on 
omF = 110*2*pi;
plot(omF,0,'*','color','r')
hold on
semilogy(omF,1e-7,'*','color','r')
semilogy(2*omF,1e-7,'*','color','r')
semilogy(3*omF,1e-7,'*','color','r')
semilogy(4*omF,1e-7,'*','color','r')
semilogy(5*omF,1e-7,'*','color','r')

%% Integrate equation of motion

%frequency content
fForc1 = 630; %350 %excite third mode (fr3 : 1750 hz)
omForc1 = 2*pi*fForc1;
TForc1 = 1/fForc1;

fForc2 = 3184; %excite third mode (fr3 : 1750 hz)
omForc2 = 2*pi*fForc2;
TForc2 = 1/fForc2;

Fext = @(t) F;%F1*cos(omForc1*t) + F2*cos(omForc2*t);

% tramp = 0.5;
% Fext = @(t) F1*((t/tramp) - (t/tramp)*heaviside(t-tramp) + heaviside(t-tramp) );

%integration time
% tint = 8*TForc1;
tint = 1/f0(1)*5;

% hnlin = TForc1/150;
hnlin = 1/f0(3)/30;

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
dynFull.nlin.disp = decodeDofsNodes(uDynFull,nNodes,2); % (node, dof of node, tsamp)
dynFull.nlin.time = TI_NL.Solution.time;

%plot solution
plot_dof_location = [1,1]; %percentage of length of beam
node2plot = find_node(plot_dof_location(1)*Lx,plot_dof_location(2)*Ly,[],nodes); % node to plot

figure
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,2,:)))


tTotAnimation = 1;
elementPlot = elements(:,1:4);
frameRate = length(dynFull.nlin.time)/tTotAnimation;
% figure
% AnimateFieldonDeformedMesh(nodes,elementPlot,uDynFull,'framerate',round(frameRate),'factor',1)

% freqContent = fft(uDynFull,nFreq,2);

%% integrate linearized model
% settings for integration
hlin = hnlin;

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
dynFull.lin.disp = decodeDofsNodes(uDynFullLin,nNodes,2); % (node, dof of node, tsamp)
dynFull.lin.time = TI_LIN.Solution.time;

plot_dof_location = [1,1]; %percentage of length of beam
node2plot = find_node(plot_dof_location(1)*Lx,plot_dof_location(2)*Ly,[],nodes); % node to plot

figure
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,2,:)))
hold on
plot(dynFull.lin.time,squeeze(dynFull.lin.disp(node2plot,2,:)))
legend('NL','LIN')

%% Construct ROMs

%put linear static solution in the basis
%uLin = Assembly.solve_system(K, Fext(0));

V1 = orth([VMs(:,1:4),SMDs]);

V2 = orth([VMs(:,1:4),PMDs1,PMDs2]);

V3 = orth([VMs(:,1:4),SMDs,PMDs1]);

Vtest = {V1,V2,V3};

% thetaMDsum = 180/pi*subspace([VMs(:,1:4),SMDs],PMDs1);
% thetaMDdiff = 180/pi*subspace([VMs(:,1:4),SMDs],PMDs2);
% 
% Z1 = orth(VMs(:,1:4));
% Z2 = orth(SMDs);
% 
% Z2Z1 = Z2'*Z1;
% [L,Sigma,R] = svd(Z2Z1);
% 
% dd = diag(Sigma);
% theta = 180/pi*acos(min(dd));

V_VMsMDs = orth([VMs(:,1:4),SMDs]);

%angle between subspaces
angle_VMsMDs_PMDsumTot = 180/pi*subspace(V_VMsMDs,PMDs1);
angle_VMsMDs_PMDdiffTot = 180/pi*subspace(V_VMsMDs,PMDs2);


angle_VMsMDs_PMDsum = zeros(size(PMDs1,2),1);
angle_VMsMDs_PMDdiff = zeros(size(PMDs2,2),1);

for ii = 1:size(PMDs1,2)
    
    angle_VMsMDs_PMDsum(ii) = 180/pi*subspace(V_VMsMDs,PMDs1(:,ii));
    angle_VMsMDs_PMDdiff(ii) = 180/pi*subspace(V_VMsMDs,PMDs2(:,ii));
    
end

thetaThreshold = 4;
maxAnglePMDsum = find(angle_VMsMDs_PMDsum > thetaThreshold);
names(maxAnglePMDsum,:)
maxAnglePMDdiff = find(angle_VMsMDs_PMDdiff > thetaThreshold);
names(maxAnglePMDdiff,:)

%% construct ROMs
V1 = orth([VMs(:,1:4),SMDs]);

V2 = orth([VMs(:,1:4),SMDs,PMDs1(:,[4 7 9]),PMDs2(:,4)]);

V3 = orth([VMs(:,1:4),PMDs2,PMDs1]);

Vtest = {V1,V2};%,V3};

%% automatic basis construction


V4 = VMs(:,1:4);
% addToV = [PMDs2,PMDs1];
addToV = [SMDs,VMs(:,5:9)];

thetaThreshold = 3;

for ii = 1:size(addToV,2)
    
    angle_VMsMDs_ii = 180/pi*subspace(V4,addToV(:,ii));
    if angle_VMsMDs_ii > thetaThreshold
         V4 = [V4,addToV(:,ii)];
         ii
    end
    
end

V4 = orth(V4);

Vtest = {V1,V4};%,V3};


%% Simulate different ROMs

dynRed = cell(length(Vtest),1);
errROM{ii} = cell(length(Vtest),1);
for ii = 1:length(Vtest)
    
    V = Vtest{ii};

    redAssembly = ReducedAssembly(myMesh,V); %reduced assembly
    redAssembly.DATA.M = V.'*Assembly.DATA.M*V;
    redAssembly.DATA.C = V.'*Assembly.DATA.C*V;  
    
    nDofRed = size(V,2);

    q0 = zeros(nDofRed,1);
    qd0 = zeros(nDofRed,1);
    qdd0 = zeros(nDofRed,1);
    
    % Instantiate object for nonlinear time integration
    TI_NL_r = ImplicitNewmark('timestep',hnlin,'alpha',0.005);
    
    % nonlinear Residual evaluation function handle
    residual_r = @(q,qd,qdd,t)residual_reduced_nonlinear(q,qd,qdd,t,redAssembly,Fext);
    
    % integrate equations with Newmark scheme
    TI_NL_r.Integrate(q0,qd0,qdd0,tint,residual_r);
    
    % project to full space
    uDynRed = V*TI_NL_r.Solution.q;
    vDynRed = V*TI_NL_r.Solution.qd;
    aDynRed = V*TI_NL_r.Solution.qdd;
    
    %save solution
    dynRed{ii}.nlin.disp = decodeDofsNodes(uDynRed,nNodes,2); % (node, dof of node, tsamp)
    dynRed{ii}.nlin.time = TI_NL_r.Solution.time;

    %compute error
    errROM{ii} = RMS_error_ROM(dynFull.nlin.disp,dynFull.nlin.time,dynRed{ii}.nlin.disp,dynRed{ii}.nlin.time,1:2,1:myMesh.nNodes);

end

%%
%displacements
plotDofLocation = [1,1]; %percentage of plate length (along x and y)
node2plot = find_node(plotDofLocation(1)*Lx,plotDofLocation(2)*Ly,[],nodes); % node to plot
dofPl = 2;

description = '';
fontsize = 15;
linewidth = 2;

figure('units','normalized','position',[0.3,0.1,0.6,0.7]); hold on
plot(dynFull.nlin.time,squeeze(dynFull.nlin.disp(node2plot,dofPl,:)),'-k','linewidth',linewidth)
plot(dynRed{1}.nlin.time,squeeze(dynRed{1}.nlin.disp(node2plot,dofPl,:)),'--r','linewidth',linewidth)
plot(dynRed{2}.nlin.time,squeeze(dynRed{2}.nlin.disp(node2plot,dofPl,:)),'-b','linewidth',linewidth)
%plot(dynRed{3}.nlin.time,squeeze(dynRed{3}.nlin.disp(node2plot,dofPl,:)),'--m','linewidth',linewidth)

xlabel('time [s]','fontsize',fontsize); ylabel('u [m]','fontsize', fontsize)
title(['normalized plotted dof location: ',num2str(plotDofLocation(1)),',',...
    num2str(plotDofLocation(2)),',  dof: ',num2str(dofPl)]);
ax = gca; grid on; ax.FontSize = fontsize; %legend('HFM','SMDs','PMDs sum diff','SMDs PMDs sum','fontsize',fontsize);
legend('HFM','SMDs','SMDs selected','PMDs all','fontsize',fontsize);%'SMDs PMDs sum','fontsize',fontsize);

% plot errror
description = '';
fontsize = 15;
linewidth = 2;

figure('units','normalized','position',[0.3,0.1,0.6,0.7]); hold on
plot(dynFull.nlin.time,errROM{1},'-k','linewidth',linewidth)
plot(dynFull.nlin.time,errROM{2},'--r','linewidth',linewidth)
%plot(dynFull.nlin.time,errROM{3},'-b','linewidth',linewidth)

xlabel('time [s]','fontsize',fontsize); ylabel('RMS error [m]','fontsize', fontsize)
title('error ');
ax = gca; grid on; ax.FontSize = fontsize; 
legend('SMDs','SMDs selected','PMDs all','fontsize',fontsize);

