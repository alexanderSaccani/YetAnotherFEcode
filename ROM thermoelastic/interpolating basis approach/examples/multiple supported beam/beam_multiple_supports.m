% author: ALEXANDER SACCANI - PHD CANDIDATE, ETH ZURICH, 
% created: 17.01.2023
% last modified: 30.01.2023

%beam with multiple supports and temperature pulse traversing the main
%beam

clearvars
clc
close all

%% Model construction

%GEOMETRY__________________________________________________________________
%main beam
l = 0.2; %length of main beam


%supports
n_supports = 2;
ls = 0.01; %lenght of supports

dl_betw_supps = l/(n_supports+1);


%SECTIONS__________________________________________________________________
%section of main beam
h = 1e-3;
b = 0.5e-2; 

%section of supports
hs = 2e-3;
bs = 2e-3;


%MATERIAL__________________________________________________________________
E       = 70e9;  % 70e9 % 200e9 % Young's modulus
rho     = 2700; % 2700 % 7850 % density
nu      = 0.3;    % nu
kappa   = 1e8; % material damping modulus 1e8
alpha_T = 23.1e-6; % thermal expansion coefficient

myThermalBeamMaterial = KirchoffMaterial();
set(myThermalBeamMaterial,'YOUNGS_MODULUS', E,'DENSITY',rho,...
    'POISSONS_RATIO',nu,'THERMAL_EXPANSION_COEFFICIENT', alpha_T);

%MESH______________________________________________________________________
%main beam
n_el_betw_supps = 15;

mbeam.n_nodes = n_el_betw_supps*(n_supports+1)+1; %number of nodes of main beam
mbeam.nelements = n_el_betw_supps*(n_supports+1); %number of element in main beam
mbeam.nodes.x = linspace(0,l,mbeam.n_nodes).'; %x nodal coordinates
mbeam.nodes.y = zeros(mbeam.n_nodes,1); %y nodal coordinates

mbeam.elconn = [(1:(mbeam.n_nodes-1)).',(2:mbeam.n_nodes).'];
mbeam.nodesID = (1:mbeam.n_nodes).';

%supports
n_el_support = 4; %number of elements for each support
n_nodes_support = n_el_support + 1;

supports = cell(n_supports,1);

for ii = 1:n_supports
    
  nodes_x_ii = ones(n_nodes_support,1)*(l/(n_supports+1))*ii;
  nodes_y_ii = flip(linspace(-ls,0,n_el_support+1)).';
  
  supports{ii}.nodes.x = nodes_x_ii;
  supports{ii}.nodes.y = nodes_y_ii;
  
  nodesID_ii = zeros(n_nodes_support,1);
  nodesID_ii(1) = n_el_betw_supps*ii + 1;
  nodesID_ii(2:end) = mbeam.n_nodes + (n_nodes_support-1)*(ii-1) + ...
      + (1:(n_nodes_support-1)).';
  
  supports{ii}.nodesID = nodesID_ii;
  
  supports{ii}.elconn = [nodesID_ii(1:end-1),nodesID_ii(2:end)]; %connectivity for elements in support
      
end

%put together the mesh of beam and supports
nNodes = mbeam.n_nodes + (n_nodes_support-1)*n_supports; 

Nodes = zeros(nNodes,2);

Nodes(1:mbeam.n_nodes,1) = mbeam.nodes.x;
Nodes(1:mbeam.n_nodes,2) = mbeam.nodes.y;

for ii = 1:n_supports
    
   Nodes(supports{ii}.nodesID(2:end),1) = supports{ii}.nodes.x(2:end);
   Nodes(supports{ii}.nodesID(2:end),2) = supports{ii}.nodes.y(2:end);
   
end

nodes_x = Nodes(:,1); %x coordinates of nodes (full structure)
nodes_y = Nodes(:,2); %y coordinates of nodes 


%create the mesh from nodes
BeamMesh = Mesh(Nodes); 

%element connectivity for supports
elconn_supports = zeros(n_el_support*n_supports,2);
for ii = 1:n_supports
    
    elconn_supports((ii-1)*n_el_support +1 : ii*n_el_support,:) = supports{ii}.elconn;
      
end

%define element constructors
ElementConstr_mbeam = @()ThermalBeamElement(b, h, myThermalBeamMaterial); %constructor for main beam
ElementConstr_supports = @()ThermalBeamElement(bs, hs, myThermalBeamMaterial); %constructor for supports

%create element table from elements
BeamMesh.create_elements_table({mbeam.elconn,elconn_supports},...
    {ElementConstr_mbeam,ElementConstr_supports});

%BOUNDARY CONDITIONS_______________________________________________________
BeamMesh.set_essential_boundary_condition([1 mbeam.n_nodes],[1 2 3],0); % main beam doubly clamped
for ii = 1:n_supports
    node_to_constr = supports{ii}.nodesID(end);
    BeamMesh.set_essential_boundary_condition(node_to_constr,[1 2 3],0);%supports are clamped
end

BeamAssembly = Assembly(BeamMesh);

%dofs free and dofs constrained
nDofsF = BeamAssembly.Mesh.EBC.nDOFs; % #dofs free (without bc)
nDofsC = nDofsF - size(BeamAssembly.Mesh.EBC.constrainedDOFs,1); % #dofs constrained (with bc)


%element connectivity matrix
nElements = mbeam.nelements + n_supports*n_el_support;

Elements = zeros(nElements,2);
Elements(1:mbeam.nelements,:) = mbeam.elconn;

for ii = 1:n_supports
    
    Elements(mbeam.nelements + (ii-1)*n_el_support +1 : mbeam.nelements ...
        + ii*n_el_support,:) = supports{ii}.elconn;
      
end

%plot mesh
figure
PlotMesh(Nodes,  Elements, 0);


%GET USEFUL QUANTITIES FROM THE CREATED MODEL______________________________
K_cold = BeamAssembly.stiffness_matrix(); % tangent stiffness matrix for T = 0, u = 0
M = BeamAssembly.mass_matrix();
C = BeamAssembly.damping_matrix();

K_C = BeamAssembly.constrain_matrix(K_cold); % should be tangent stiffness matrix for T = 0, u = 0
M_C = BeamAssembly.constrain_matrix(M);

%% Modal Analysis of the cold structure

n_VMs = 5; % first n_VMs modes with lowest frequency calculated 

[VM,omega2] = eigs(K_C,M_C,n_VMs,'SM'); 
omega = sqrt(diag(omega2));
%normalize eigenvectors
for ii = 1:n_VMs
    VM(:,ii) = VM(:,ii)./vecnorm(VM(:,ii));
end
VM = BeamAssembly.unconstrain_vector(VM); %vibration modes

%decode output
nDofPerNode = 3;
VM = decodeDofsNodes(VM,nNodes,nDofPerNode); %(#node,#dof,#VM)

%plot VMs
mode2plot = [1,2,3];

for ii = 1:length(mode2plot)

jj = mode2plot(ii);

figure
PlotMesh(Nodes,  Elements, 0); hold on;
PlotFieldonDeformedMesh(Nodes,Elements,[VM(:,1,jj),VM(:,2,jj)],'factor',1, 'color','r');

title(['VM ',num2str(mode2plot(ii)),', Frequency = '...
     num2str(omega(mode2plot(ii))/(2*pi)) ' Hz'])

end

% %plot VMs
% mode2plot = [1,2];
% for ii = 1:length(mode2plot)
%   
%     VM_ii = squeeze(VM(:,:,mode2plot(ii)));
%     
%     u_ii = VM_ii(:,1); % long. displ.
%     v_ii = VM_ii(:,2); % transv. displ
%     r_ii = VM_ii(:,3); % rotation
%     
%     figure
%     subplot 411
%     plot(nodes_x,nodes_y,'.-','markersize',10); hold on;
%     plot(nodes_x + u_ii, nodes_y + v_ii,'.-','markersize',10)
%     title(['VM ',num2str(mode2plot(ii)),', Frequency = '...
%     num2str(omega(mode2plot(ii))/(2*pi)) ' Hz'])
%     
%     subplot 412
%     plot(nodes_x,v_ii,'-'); title('tranverse disp.'); xlabel('x')
%     subplot 413
%     plot(nodes_x,u_ii,'-'); title('longitudinal disp.'); xlabel('x')
%     subplot 414
%     plot(nodes_x,r_ii,'-'); title('rotation'); xlabel('x')
%     
% end

%% Define Temperature distribution profile (shape of profile)

%define T distribution (shape)
T_ampl = 100; %amplitude of thermal pulse
k = 1; %factor that defines, with the number of supports, the width of the pulse
p = k*l/(n_supports+1); %width of sine pulse
x0 = @(xc) xc - p/2; %left extreme of pulse
T_dyn_p = @(xc) T_ampl*(sin(pi*(nodes_x-x0(xc))/p)).^2.*(heaviside(nodes_x-x0(xc)) - ...
    heaviside((nodes_x-x0(xc))-p)); %define the temperature as a function of xc only (T profile parametrized with xc (centre of pulse))


%% Modal analysis of the hot structure

%define sampling (samples created by varying xc)
ecc = -0.5; %eccentricity of applied pulse center. fraction of distance between supports w.r.t. the center of the pulse is shifted
dl_ecc = ecc*dl_betw_supps;
p_sampl = [-p/2,dl_betw_supps+dl_ecc:dl_betw_supps:l-1e-15*l]; %compute VMs for cold beam, and for other T configs

%generate corresponding temperature profiles
T_sampls = zeros(length(nodes_x),length(p_sampl),1);
for ii = 1:length(p_sampl)
    T_sampls(:,ii) = T_dyn_p(p_sampl(ii));
end

VMs_at_eq = 1; %do you want to compute VMs at thermal equilibrium? (set to 1 if yes, otherwise set to 0)
how_many_VMs = 4;
th_modes_analysis = VMMD_thermal(BeamAssembly,T_sampls,VMs_at_eq,how_many_VMs,1);

VMs_hot = th_modes_analysis.VMs;
static_eq = th_modes_analysis.eq;
omega = th_modes_analysis.omega;
MDs_hot = th_modes_analysis.MDs;

color_list = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE',...
    '#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};

mode2plot = [1,2];
T_sampl_2_plot = 1:1:length(p_sampl); %plot all T distributions

%VMs
scale_factor = 2;
for ii = 1:length(mode2plot)
    
    VM2pl = mode2plot(ii);
    figure('units','normalized','position',[.3 .3 .4 .4]); hold on;
    title(['VM ',num2str(mode2plot(ii))]);
    
    for jj = 1:length(T_sampl_2_plot)
        %decode output
        VM_iijj = decodeDofsNodes(VMs_hot{jj}(:,VM2pl),nNodes,nDofPerNode);
        PlotFieldonDeformedMesh(Nodes,Elements,[VM_iijj(:,1),VM_iijj(:,2)],'factor',scale_factor, 'color',color_list{jj});
        PlotFieldonDeformedMesh(Nodes,Elements,[-VM_iijj(:,1),-VM_iijj(:,2)],'factor',scale_factor, 'color',color_list{jj});
    end

end

%Tsamples
figure('units','normalized','position',[.3 .3 .4 .4]); hold on;
title('T distribution samples'); xlabel('x [m]'); ylabel('T [K]');
for jj = 1:length(T_sampl_2_plot)
    plot( nodes_x(1:mbeam.n_nodes),T_sampls(1:mbeam.n_nodes,T_sampl_2_plot(jj)),'color',color_list{jj})
end

%static equilibrium
static_eq = decodeDofsNodes(static_eq,nNodes,nDofPerNode);
scale_factor = 10;
figure('units','normalized','position',[.3 .3 .4 .4]); hold on;
for jj = 1:length(T_sampl_2_plot)
   title(['Static thermal equilibrium, scale factor: ',num2str(scale_factor)]); 
   PlotFieldonDeformedMesh(Nodes,Elements,[static_eq(:,1,jj),static_eq(:,2,jj)],'factor',scale_factor, 'color',color_list{jj});
end

%MDs
scale_factor = 2;
for ii = 1:size(MDs_hot{1,1},2)
    
    MD2pl = ii;
    figure('units','normalized','position',[.3 .3 .4 .4]); hold on;
    title(['MD ',num2str(ii)]);
    
    for jj = 1:length(T_sampl_2_plot)
        %decode output
        MD_iijj = decodeDofsNodes(MDs_hot{jj}(:,MD2pl),nNodes,nDofPerNode);
        PlotFieldonDeformedMesh(Nodes,Elements,[MD_iijj(:,1),MD_iijj(:,2)],'factor',scale_factor, 'color',color_list{jj});
        PlotFieldonDeformedMesh(Nodes,Elements,[-MD_iijj(:,1),-MD_iijj(:,2)],'factor',scale_factor, 'color',color_list{jj});
    end

end


%display natural frequencies values
disp('natural frequencies are: ');
disp(omega)


%% Static analysis

%DEFINE LOAD_______________________________________________________________
%concentrated load
Fc_mag = -70; %magnitude of load 
Fc_n_supp = 1; %define after which support to put the concentrated force (must be less or equal to the number of supports) 
%(Fc_n_supp must be less or equal to n_supports)
if Fc_n_supp > n_supports
    error('incorrect input')
end
Fc_ecc = 0.3; %define where to put the load (range is from 0 to 1)
Fc_dof = 2; %which dof to excite

Fc_loc = Fc_n_supp*dl_betw_supps + Fc_ecc*dl_betw_supps;
fext_node = find_node(Fc_loc,0,[],[mbeam.nodes.x,mbeam.nodes.y]); % node where to put the force
fext_dofs = get_index(fext_node, nDofPerNode);

Fc = zeros(nDofsF,1); %concentrated force
Fc(fext_dofs(Fc_dof)) = Fc_mag; 

F = Fc;

%DEFINE CENTER OF T DISTRIBUTION___________________________________________
ecc = +0.5; %eccentricity of applied pulse center. fraction of distance between supports w.r.t. the center of the pulse is shifted
dl_ecc = ecc*dl_betw_supps;

xc_n_supp = 1; %define after which support to put the center of T distr (must be less or equal to the number of supports) 
if xc_n_supp > n_supports
    error('incorrect input')
end
xc_st = xc_n_supp*dl_betw_supps + dl_ecc;
T_st =  T_dyn_p(xc_st); %static nodal T profile

%RUN LINEAR AND NONLINEAR STATIC ANALYSIS__________________________________
% linear static solution
[K_hot,gth] = BeamAssembly.tangent_stiffness_and_force(zeros(nDofsF,1),T_st);
u_lin = BeamAssembly.solve_system(K_hot, F-gth);   %subtract the internal force generated  by the temperature
static.lin = reshape(u_lin,3,[]).';	% Linear response

% nonlinear static solution
BeamAssembly.DATA.K = K_hot;    % this value of K is used only for first guess solution 
[~,u] = static_equilibrium_thermal(BeamAssembly, F, T_st,'display', 'iter-detailed');
static.nlin = reshape(u,nDofPerNode,[]).';     % Nonlinear response


%PLOT RESULTS OF STATIC ANALYSIS___________________________________________
scale_factor = 1;

figure('units','normalized','position',[.1 .1 .8 .8])
subplot 211; 
plot( nodes_x(1:mbeam.n_nodes),T_st(1:mbeam.n_nodes))
xlabel('x [m]'); ylabel('T [K]');
title('Static analysis. T distribution')
subplot 212;
PlotMesh(Nodes,  Elements, 0); hold on;
PlotFieldonDeformedMesh(Nodes,Elements,[static.lin(:,1),static.lin(:,2)],'factor',scale_factor, 'color','b');
PlotFieldonDeformedMesh(Nodes,Elements,[static.nlin(:,1),static.nlin(:,2)],'factor',scale_factor, 'color','r');
legend('undeformed','linear','nonlinear')
title('Static analysis, deformed shape')

%% Dynamic Analysis
%DEFINE LOAD_______________________________________________________________
%frequncy of load
om_fext = (omega(1)+omega(2))/2; 
T_fext = 2*pi/om_fext; % period of ex. forcing

%shape of load
Fc_mag = -30; %magnitude of load 
Fc_n_supp = 0; %define after which support to put the concentrated force (must be less or equal to the number of supports) 
%(Fc_n_supp must be less or equal to n_supports)
if Fc_n_supp > n_supports
    error('incorrect input')
end
Fc_ecc = 0.3; %define where to put the load (range is from 0 to 1)
Fc_dof = 2; %which dof to excite

%construct load
Fc_loc = Fc_n_supp*dl_betw_supps + Fc_ecc*dl_betw_supps;
fext_node = find_node(Fc_loc,0,[],[mbeam.nodes.x,mbeam.nodes.y]); % node where to put the force
fext_dofs = get_index(fext_node, nDofPerNode);

Fc = zeros(nDofsF,1); %concentrated force
Fc(fext_dofs(Fc_dof)) = Fc_mag; 

F = Fc;

F_ext = @(t) F*sin(om_fext*t); %external forcing

%DEFINE TEMPERATURE VARIATION______________________________________________
%define how center of pulse varies in time
A = l+p;  %amplitude of oscillation of pulse
xci = -p/2; %center location at initial time
eps = 0.001; %om_fex/om_th (th = thermal forcing)

%dependent parameters
om_th = eps*om_fext;  % angular frequency of thermal mode oscillation in time
T_th = 2*pi/om_th; % period of temp. oscillations

xc_dyn = @(t) xci + A*sin(om_th*t); %define how center of pulse xc variates in time


T_dyn_t = @(t) T_dyn_p(xc_dyn(t)); %evaluate the temperature profile in time


% %plot temperature profile (animation)
% time = linspace(0,T_th/4,100);
% figure
% for ii = 1:length(time)
%     plot(nodes_x,T_dyn_t(time(ii)));
%     axis([0,l,0,T_ampl]);
%     pause(0.03);
% end

%DEFINE INTEGRATION INTERVAL (the same for all models)_____________________
t_int = T_th/4; % integration interval 

%INTEGRATE LINEARIZED MODEL________________________________________________
%settings for integration
h = T_fext/50; % time step for integration

%Precompute data for Assembly object
BeamAssembly.DATA.M = M;
BeamAssembly.DATA.K = K_cold;
BeamAssembly.DATA.C = C;

%Initial condition: equilibrium
q0 = zeros(nDofsC,1);
qd0 = zeros(nDofsC,1);
qdd0 = zeros(nDofsC,1);

%Instantiate object for linear time integration
TI_LIN = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

%nonlinear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_thermal_linear(q,qd,qdd,t,...
    BeamAssembly, F_ext, T_dyn_t);

%integrate equations with Newmark scheme
TI_LIN.Integrate(q0,qd0,qdd0,t_int,residual);

%store solution in more suitable array (only displacements)
u_dyn = BeamAssembly.unconstrain_vector(TI_LIN.Solution.q); 
dyn.lin.disp = decodeDofsNodes(u_dyn,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn.lin.time = TI_LIN.Solution.time;

%INTEGRATE NONLINEAR MODEL_________________________________________________
% settings for integration
h = T_fext/50; % time step for integration

%Precompute data for Assembly object
BeamAssembly.DATA.M = M;
BeamAssembly.DATA.K = K_cold;
BeamAssembly.DATA.C = C;

%Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

%nonlinear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_thermal_nonlinear(q,qd,qdd,t,...
    BeamAssembly,F_ext, T_dyn_t);

%integrate equations with Newmark scheme
TI_NL.Integrate(q0,qd0,qdd0,t_int,residual);

%store solution in more suitable array (only displacements)
u_dyn = BeamAssembly.unconstrain_vector(TI_NL.Solution.q); 
dyn.nlin.disp = decodeDofsNodes(u_dyn,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn.nlin.time = TI_NL.Solution.time;

%QUASISTATIC SOLUTION______________________________________________________
dyn.quasistatic.time = linspace(0,t_int,1e2); %define time vector
u_quasist = quasistatic_solution( BeamAssembly, T_dyn_t, dyn.quasistatic.time ); %compute quasistatic solution (F=0)
dyn.quasistatic.disp = decodeDofsNodes(u_quasist,nNodes,nDofPerNode); % (node, dof of node, tsamp)

%%
%PLOT RESULTS OF DYNAMIC RESPONSE__________________________________________
plot_dof_loc = 0.5; %percentage of length of beam
node2plot = find_node(plot_dof_loc*l,0,[],[mbeam.nodes.x,mbeam.nodes.y]); % node to plot
%node2plot = fext_node; % uncomment to plot time history of forced node

figure('units','normalized','position',[.1 .1 .8 .8]);
subplot 311; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,1,:)),'-','color','r','linewidth',2);
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,1,:)),'--','color','cyan','linewidth',2);
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,1,:)),'--','color','green','linewidth',2);
xlabel('t [s]'); ylabel('u [m]'); grid on; axis tight;  legend('nonlinear','linear','quasistatic','linewidth',2);
title(['axial displacements (u), transversal displacements (v) and rotations (w) at ', num2str(plot_dof_loc),'l'])

subplot 312; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,2,:)),'-','color','r','linewidth',2);
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,2,:)),'--','color','cyan','linewidth',2);
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,2,:)),'--','color','green','linewidth',2);
xlabel('t [s]'); ylabel('v [m]'); grid on; axis tight; legend('nonlinear','linear','quasistatic','linewidth',2);

subplot 313; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,3,:)),'-','color','r','linewidth',2);
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,3,:)),'--','color','cyan','linewidth',2);
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,3,:)),'--','color','green','linewidth',2);
xlabel('t [s]'); ylabel('w [-]'); grid on; axis tight;  legend('nonlinear','linear','quasistatic','linewidth',2);

%% Reduced Order Model Construction (varying basis)
% define number of equidistant samples
n_ROMs = 20;

xc_sampls = linspace(-p/2,l+p/2,n_ROMs)'; %sample locations of center thermal pulse

T_sampls = zeros(nNodes,length(xc_sampls)); % T nodal distribution samples

%generate corresponding temperature profiles
for ii = 1:length(xc_sampls)
    T_sampls(:,ii) = T_dyn_p(xc_sampls(ii));
end

%construct model___________________________________________________________
% define how many VMs to use in modal reduction (also the corresponding
% modal derivatives are added)
number_VMs = 6;

ROMs.fullAssembly = BeamAssembly;
ROMs.parameters = xc_sampls'; 
ROMs.models = multiple_ROMs_thermal(BeamAssembly, T_sampls, number_VMs); 
% cell array, with in each cell a Model. Each model is a struct array with field V (for now, then could be extended to M, K, K2, K3)

%% Integrate nonlinear ROM (Varying basis)
% settings for integration
h = T_fext/50; % time step for integration

% Initial condition: equilibrium
nDofRed = size(ROMs.models{1}.V,2);

q0 = zeros(nDofRed,1);
qd0 = zeros(nDofRed,1);
qdd0 = zeros(nDofRed,1);

% Instantiate object for nonlinear time integration
TI_NL_r = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear Residual evaluation function handle
residual_r = @(q,qd,qdd,t)residual_ROM_thermal(q,qd,qdd,t,ROMs,F_ext,xc_dyn,T_dyn_p,u_quasist,dyn.quasistatic.time);

% integrate equations with Newmark scheme
TI_NL_r.Integrate(q0,qd0,qdd0,t_int-1e-3*t_int,residual_r);

% postprocess output 
u_dyn_r = reduced_to_full_thermal(TI_NL_r.Solution.q,TI_NL_r.Solution.time,ROMs,xc_dyn,u_quasist,dyn.quasistatic.time);

u_dyn_r = BeamAssembly.unconstrain_vector(u_dyn_r); 
dyn_r.nlin.disp = decodeDofsNodes(u_dyn_r,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn_r.nlin.time = TI_NL_r.Solution.time;

%quasistatic response
u_qs = quasistatic_from_ROMs(dyn.nlin.time,ROMs,xc_dyn);
u_qs = BeamAssembly.unconstrain_vector(u_qs);
u_qs = decodeDofsNodes(u_qs,nNodes,nDofPerNode);

%% Integrate nonlinear ROM (constant basis)
% The ROM with constant basis, can be obtained (neglecting computational
% efficiency) by setting all the basis of the different precomputed ROMs to
% the basis of the cold structure (the first one). In this way the rest of
% the code remains the same.    N.B. NOT GOOD TO CHECK EFFICIENCY!
ROMs_cbasis = ROMs;
V_const = ROMs.models{1}.V;
for ii = 1:n_ROMs
    ROMs_cbasis.models{ii}.V = V_const;
end

% settings for integration
h = T_fext/50; % time step for integration

% Initial condition: equilibrium
nDofRed = size(ROMs.models{1}.V,2);

q0 = zeros(nDofRed,1);
qd0 = zeros(nDofRed,1);
qdd0 = zeros(nDofRed,1);

% Instantiate object for nonlinear time integration
TI_NL_r_cbasis = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear Residual evaluation function handle
residual_r_cbasis = @(q,qd,qdd,t)residual_ROM_thermal(q,qd,qdd,t,ROMs_cbasis,F_ext,xc_dyn,T_dyn_p,u_quasist,dyn.quasistatic.time);

% integrate equations with Newmark scheme
TI_NL_r_cbasis.Integrate(q0,qd0,qdd0,t_int-1e-3*t_int,residual_r_cbasis);

% postprocess output 
u_dyn_r_cbasis = reduced_to_full_thermal(TI_NL_r_cbasis.Solution.q,TI_NL_r_cbasis.Solution.time,ROMs_cbasis,xc_dyn,u_quasist,dyn.quasistatic.time);

u_dyn_r_cbasis = BeamAssembly.unconstrain_vector(u_dyn_r_cbasis); 
dyn_r_cb.nlin.disp = decodeDofsNodes(u_dyn_r_cbasis,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn_r_cb.nlin.time = TI_NL_r_cbasis.Solution.time;

%% Plot comparison
% plot results of nonlinear dynamic analysis_______________________________
plot_dof_loc = 0.3; %percentage of length of beam
node2plot = find_node(plot_dof_loc*l,0,[],[mbeam.nodes.x,mbeam.nodes.y]); % node to plot

figure('units','normalized','position',[.1 .1 .8 .8]); 
subplot 311; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,1,:)),'-','color','r','linewidth',2);
plot(dyn_r.nlin.time,squeeze(dyn_r.nlin.disp(node2plot,1,:)),'--','color','b','linewidth',2);
plot(dyn_r_cb.nlin.time,squeeze(dyn_r_cb.nlin.disp(node2plot,1,:)),'-.','color','m','linewidth',2);
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,1,:)),'--','color','cyan','linewidth',2);
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,1,:)),'--','color','g','linewidth',2);
plot(dyn.nlin.time,squeeze(u_qs(node2plot,1,:)),'-.','color','k','linewidth',2)
xlabel('t [s]'); ylabel('u [m]'); grid on; axis tight; legend('full','ROM','ROM const V','lin','quasistatic','expansion point');
title(['comparison of ROM and HFM solutions, dof at ',num2str(plot_dof_loc*100), '% of beam length'])
axis([0,t_int,min(min(dyn_r_cb.nlin.disp(node2plot,1,:))),max(max(dyn_r_cb.nlin.disp(node2plot,1,:)))]);

subplot 312; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,2,:)),'-','color','r','linewidth',2);
plot(dyn_r.nlin.time,squeeze(dyn_r.nlin.disp(node2plot,2,:)),'--','color','b','linewidth',2);
plot(dyn_r_cb.nlin.time,squeeze(dyn_r_cb.nlin.disp(node2plot,2,:)),'-.','color','m','linewidth',2);
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,2,:)),'--','color','cyan','linewidth',2);
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,2,:)),'--','color','g','linewidth',2);
plot(dyn.nlin.time,squeeze(u_qs(node2plot,2,:)),'-.','color','k','linewidth',2)
xlabel('t [s]'); ylabel('v [m]'); grid on; axis tight; legend('full','ROM','ROM const V','lin','quasistatic','expansion point');
axis([0,t_int,min(min(dyn_r_cb.nlin.disp(node2plot,2,:))),max(max(dyn_r_cb.nlin.disp(node2plot,2,:)))]);

subplot 313; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,3,:)),'-','color','r','linewidth',2);
plot(dyn_r.nlin.time,squeeze(dyn_r.nlin.disp(node2plot,3,:)),'--','color','b','linewidth',2);
plot(dyn_r_cb.nlin.time,squeeze(dyn_r_cb.nlin.disp(node2plot,3,:)),'-.','color','m','linewidth',2);
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,3,:)),'--','color','cyan','linewidth',2);
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,3,:)),'--','color','g','linewidth',2);
plot(dyn.nlin.time,squeeze(u_qs(node2plot,3,:)),'-.','color','k','linewidth',2)
xlabel('t [s]'); ylabel('w [-]'); grid on; axis tight; legend('full','ROM','ROM const V','lin','quasistatic','expansion point');
axis([0,t_int,min(min(dyn_r_cb.nlin.disp(node2plot,3,:))),max(max(dyn_r_cb.nlin.disp(node2plot,3,:)))]);





















































