% ALEXANDER SACCANI - PHD CANDIDATE, ETH ZURICH, 
% created: 1.12.2022
% last modified: 1.12.2022

% Straight temperature dependent beam template


% single campled, axial loading

clearvars ; close all; clc

%% Model construction

%GEOMETRY__________________________________________________________________
l = 0.1; 
h = 1e-3;
b = 1e-2; 
w = 0.1; % this parameter varies from 0 to 0.5 (straight beam to half circumference)


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
nElements = 60;

%define nodes
R = (l/4 + l*w^2)/(2*w);
yy = R - l*w;
xx = l/2;
th0 = atan(yy/xx);

th = linspace(th0,pi-th0,nElements+1);

nodes_x = (R*cos(th) + l/2)'; %x coordinates for nodes
nodes_y = (R*sin(th) - (R-l*w))'; %y coordinates for nodes

nNodes = length(nodes_x); %number of nodes
nDofPerNode = 3;

%plot nodes
figure; title('beam mesh')
plot(nodes_x,nodes_y,'*-');
axis('equal')

%create mesh from nodes
Nodes = [nodes_x, nodes_y];
BeamMesh = Mesh(Nodes); 

%define elements (connectivity matrix)
Elements = [1:nNodes-1;2:nNodes].'; % ROW -> # element, 
                                    % COLUMN -> # of node connected
myElementConstructor = @()ThermalBeamElement(b, h, myThermalBeamMaterial); %type of element

%create element table from elements
BeamMesh.create_elements_table(Elements,myElementConstructor);


%BOUNDARY CONDITIONS_______________________________________________________
BeamMesh.set_essential_boundary_condition([1 BeamMesh.nNodes],[1 2 3],0) % Doubly clamped beam
%BeamMesh.set_essential_boundary_condition(1,[1 2 3],0) % clamped beam  (node,[constrained dofs],imposed disp.)

BeamAssembly = Assembly(BeamMesh);

nDofsF = BeamAssembly.Mesh.EBC.nDOFs; % #dofs free (without bc)
nDofsC = nDofsF - size(BeamAssembly.Mesh.EBC.constrainedDOFs,1); % #dofs constrained (with bc)


%GET USEFUL QUANTITIES FROM THE CREATED MODEL______________________________
K_cold = BeamAssembly.stiffness_matrix(); %tangent stiffness matrix for T = 0, u = 0
M = BeamAssembly.mass_matrix();
C = BeamAssembly.damping_matrix();

K_cold_C = BeamAssembly.constrain_matrix(K_cold); % should be tangent stiffness matrix for T = 0, u = 0
M_C = BeamAssembly.constrain_matrix(M);

%% Modal Analysis of the cold structure

n_VMs = 2; % first n_VMs modes with lowest frequency calculated 

[VM,omega2] = eigs(K_cold_C,M_C,n_VMs,'SM'); 
omega = sqrt(diag(omega2));
%normalize with respect to max. amplitude
for ii = 1:n_VMs
    VM(:,ii) = VM(:,ii)/max(sqrt(sum(VM(:,ii).^2,2)));
end
VM = BeamAssembly.unconstrain_vector(VM); %vibration modes

%decode output
VM = decodeDofsNodes(VM,nNodes,nDofPerNode);

%plot VMs
mode2plot = [1,2];
for ii = 1:length(mode2plot)
  
    VM_ii = squeeze(VM(:,:,mode2plot(ii)));
    
    u_ii = VM_ii(:,1); % long. displ.
    v_ii = VM_ii(:,2); % transv. displ
    r_ii = VM_ii(:,3); % rotation
    
    figure
    subplot 411
    plot(nodes_x,nodes_y,'.-','markersize',10); hold on;
    plot(nodes_x + u_ii, nodes_y + v_ii,'.-','markersize',10)
    title(['VM ',num2str(mode2plot(ii)),', Frequency = '...
    num2str(omega(mode2plot(ii))/(2*pi)) ' Hz  (cold)'])
    
    subplot 412
    plot(nodes_x,v_ii,'-'); title('tranverse disp.'); xlabel('x')
    subplot 413
    plot(nodes_x,u_ii,'-'); title('longitudinal disp.'); xlabel('x')
    subplot 414 
    plot(nodes_x,r_ii,'-'); title('rotation'); xlabel('x')
    
end


%% Define Loads
% Nodal force
F_c = zeros(nDofsF,1);  %concentrated load

% transverse load at midspan
fext_node = find_node(l/2,0,[],Nodes); % node where to put the force
fext_dofs = get_index(fext_node, nDofPerNode);
F_c(fext_dofs(2)) = 3e2;%3e2;%4e3; %5e3

%uniform body force
Pressure_v = 0;%1e3;%+3e3; %0.5e4;%1e2; % load per unit length in the transversal direction
F_v_uniform = Pressure_v*BeamAssembly.uniform_body_force();
F_v_uniform(1:3:end) = 0*F_v_uniform(1:3:end); %set to 0 the forces 
% corresponding to the horizontal dofs (pressure acts on the vertical direction only)

%define total external load
F = F_c + F_v_uniform;

% % Nodal force
% % % F = zeros(nDofsF,1);
% % % 
% % % % transverse load at midspan
% % % nodeF = find_node(l/2,0,[],Nodes); % node where to put the force
% % % node_force_dofs = get_index(nodeF, nDofPerNode );
% % % F(node_force_dofs(2)) = 1; 

%% Static Analysis 
% Define static thermal sin^2 distribution
T_ampl = 100;
p = l/10; %width of thermal pulse

xc_st = l/2; %center of thermal pulse
x0_st = xc_st - p/2; %left extreme of pulse
T_st = T_ampl*(sin(pi*(nodes_x-x0_st)/p)).^2.*(heaviside(nodes_x-x0_st) - ...
    heaviside((nodes_x-x0_st)-p)); %static temp profile (pulse centered at l/2)

% linear static solution___________________________________________________
[K_hot,gth] = BeamAssembly.tangent_stiffness_and_force(zeros(nDofsF,1),T_st);
u_lin = BeamAssembly.solve_system(K_hot, F-gth);   %subtract the internal force generated  by the temperature
static.lin = reshape(u_lin,3,[]).';	% Linear response

% nonlinear static solution (here instead the temperature is considered)
BeamAssembly.DATA.K = K_hot;    % this value of K is used only for first guess solution (T not considered in the linear guess!)
[~,u] = static_equilibrium_thermal(BeamAssembly, F, T_st,'display', 'iter-detailed');
static.nlin = reshape(u,nDofPerNode,[]).';     % Nonlinear response

%plot static displacements_________________________________________________
scl = 1;   %scaling factor

figure('units','normalized','position',[.1 .1 .8 .8])
subplot 411; 
plot(nodes_x,T_st); xlabel('x [m]'); ylabel('T [K]'); 
title('Static analysis. T distribution')

subplot 412; 
plot(nodes_x,static.nlin(:,1),'color','r'); hold on;
plot(nodes_x,static.lin(:,1),'--', 'color','cyan')
xlabel('x [m]'); ylabel('u [m]'); title('long. displ');legend('nonlinear','linear')

subplot 413; 
plot(nodes_x,static.nlin(:,2),'color','r'); hold on;
plot(nodes_x,static.lin(:,2),'--', 'color','cyan');
xlabel('x [m]'); ylabel('v [m]'); title('tranv. displ'); legend('nonlinear','linear')

subplot 414; 
plot(nodes_x,nodes_y,'.-','markersize',10,'color','b'); hold on;
plot(nodes_x + scl*static.nlin(:,1), nodes_y + scl*static.nlin(:,2), '.-', 'markersize',10,'color','r')
plot(nodes_x + scl*static.lin(:,1), nodes_y + scl*static.lin(:,2), '--', 'color','cyan')
title('static deformation shape'); xlabel('x [m]'); legend('undeformed','nonlinear','linear')
axis('equal')

% %plot static displacements
% ud = static.nlin(:,1);
% vd = static.nlin(:,2);
% wd = static.nlin(:,3);
% scl = 50;   %scaling factor
% 
% figure('units','normalized','position',[.1 .1 .8 .8])
% subplot 411; plot(nodes_x,T_st); xlabel('x'); ylabel('T'); 
% title('Static analysis. T distribution')
% subplot 412; plot(nodes_x,ud); xlabel('x'); ylabel('u'); title('long. displ');
% subplot 413; plot(nodes_x,vd); xlabel('x'); ylabel('v'); title('tranv. displ');
% subplot 414; plot(nodes_x,nodes_y,'.-','markersize',10,'color','b'); hold on;
% plot(nodes_x + scl*ud, nodes_y + scl*vd, '.-', 'markersize',10,'color','r')
% title('static deformation shape'), xlabel('x')


%% Dynamic response using Implicit Newmark

%force and thermal load____________________________________________________

%define external forcing
om_fext = (omega(1)+omega(2))/2;%(omega(1) + omega(2))/2;% + (omega(2) - omega(1))/12;
T_fext = 2*pi/om_fext; % period of ex. forcing

%construct dynamic forcing vector
F_ext = @(t) F*sin(om_fext*t);

%define thermal load 
T_ampl = 100;
p = l/10; %width of thermal pulse
A = l+p;  %amplitude of oscillation of pulse
xci = -p/2; %center location at initial time
eps = 0.01; %om_fex/om_th (th = thermal forcing)
om_th = eps*om_fext;  % angular frequency of thermal mode oscillation in time

T_th = 2*pi/om_th; % period of temp. oscillations

x0 = @(xc) xc - p/2; %left extreme of pulse
T_dyn_xc = @(xc) T_ampl*(sin(pi*(nodes_x-x0(xc))/p)).^2.*(heaviside(nodes_x-x0(xc)) - ...
    heaviside((nodes_x-x0(xc))-p)); %define the temperature as a function of xc only (T profile parametrized with xc)

xc_dyn = @(t) xci + A*sin(om_th*t); %define how center of pulse xc variates in time
T_dyn_t = @(t) T_dyn_xc(xc_dyn(t)); %evaluate the temperature profile in time

% %plot temperature profile
% time = linspace(0,T_th/4,100);
% figure
% for ii = 1:length(time)
%     plot(nodes_x,T_dyn_t(time(ii)));
%     axis([0,l,0,T_ampl]);
%     pause(0.03);
% end

% Integrate linearized model_______________________________________________
% settings for integration
tmax = T_th/4; % integration interval
h = T_fext/50; % time step for integration

% Precompute data for Assembly object
BeamAssembly.DATA.M = M;
BeamAssembly.DATA.K = K_cold;
BeamAssembly.DATA.C = C;

% Initial condition: equilibrium
q0 = zeros(nDofsC,1);
qd0 = zeros(nDofsC,1);
qdd0 = zeros(nDofsC,1);

% Instantiate object for linear time integration
TI_LIN = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% nonlinear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_thermal_linear(q,qd,qdd,t,...
    BeamAssembly, F_ext, T_dyn_t);

% integrate equations with Newmark scheme
TI_LIN.Integrate(q0,qd0,qdd0,tmax,residual);

% store solution in more suitable array (only displacements)
u_dyn = BeamAssembly.unconstrain_vector(TI_LIN.Solution.q); 
dyn.lin.disp = decodeDofsNodes(u_dyn,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn.lin.time = TI_LIN.Solution.time;

% Integrate nonlinear EOM__________________________________________________

% settings for integration
tmax = T_th/4; % integration interval
h = T_fext/50; % time step for integration

% Precompute data for Assembly object
BeamAssembly.DATA.M = M;
BeamAssembly.DATA.K = K_cold;
BeamAssembly.DATA.C = C;

% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_thermal_nonlinear(q,qd,qdd,t,...
    BeamAssembly,F_ext, T_dyn_t);

% integrate equations with Newmark scheme
TI_NL.Integrate(q0,qd0,qdd0,tmax,residual);

% store solution in more suitable array (only displacements)
u_dyn = BeamAssembly.unconstrain_vector(TI_NL.Solution.q); 
dyn.nlin.disp = decodeDofsNodes(u_dyn,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn.nlin.time = TI_NL.Solution.time;

% Quasistatic Solution_____________________________________________________
dyn.quasistatic.time = linspace(0,tmax,1e2); %define time vector
u_quasist = quasistatic_solution( BeamAssembly, T_dyn_t, dyn.quasistatic.time ); %compute quasistatic solution (F=0)
dyn.quasistatic.disp = decodeDofsNodes(u_quasist,nNodes,nDofPerNode); % (node, dof of node, tsamp)

% plot results of nonlinear dynamic analysis_______________________________
plot_dof_location = 0.3; %percentage of length of beam
node2plot = find_node(plot_dof_location*l,0,[],Nodes); % node where to put the force

figure('units','normalized','position',[.1 .1 .8 .8]);
subplot 311; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,1,:)),'-','color','r');
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,1,:)),'--','color','cyan');
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,1,:)),'--','color','green');
xlabel('t [s]'); ylabel('u [m]'); grid on; axis tight;  legend('nonlinear','linear','quasistatic');

subplot 312; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,2,:)),'-','color','r');
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,2,:)),'--','color','cyan');
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,2,:)),'--','color','green');
xlabel('t [s]'); ylabel('v [m]'); grid on; axis tight; legend('nonlinear','linear','quasistatic');

subplot 313; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,3,:)),'-','color','r');
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,3,:)),'--','color','cyan');
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,3,:)),'--','color','green');
xlabel('t [s]'); ylabel('w [-]'); grid on; axis tight;  legend('nonlinear','linear','quasistatic');


% %% plot results in time (animation)______________________________________
% time_tmp = dyn.nlin.time; % time vec.
% u_tmp = squeeze(dyn.nlin.disp(:,1,:)); % long. displ.
% v_tmp = squeeze(dyn.nlin.disp(:,2,:)); % tranv. displ.
% 
% n_samp_plot = 500;
% dind_plot = round(length(time_tmp)/n_samp_plot);
% ind_plot = 1:dind_plot:length(time_tmp);
% time_tmp = time_tmp(ind_plot); % time vec.
% u_tmp = u_tmp(:,ind_plot); % long. displ.
% v_tmp = v_tmp(:,ind_plot);
% 
% u_ax_bounds = [0,l,1.1*min(min(u_tmp)),1.1*max(max(u_tmp))];
% v_ax_bounds = [0,l,1.1*min(min((v_tmp))),1.1*max(max(v_tmp))];
% 
% figure('units','normalized','position',[.1 .1 .8 .8]);
% for ii = 1:length(time_tmp)
%    
%    time_ii = time_tmp(ii);
%    
%    subplot 311
%    plot(nodes_x,T_dyn(time_tmp(ii)));
%    axis([0,l,0,T_ampl]);
%    xlabel('x'); ylabel('T');
%    
%    subplot 312
%    plot(nodes_x, u_tmp(:,ii))
%    axis(u_ax_bounds)
%    xlabel('x'); ylabel('u')
%    
%    subplot 313
%    plot(nodes_x, v_tmp(:,ii))
% %  axis([0,l,-1,1])
%    axis(v_ax_bounds)
%    xlabel('x'); ylabel('v')
%    
%    pause(0.001);
%    
% end

%% Thermal VMs analysis (just to visualize thermal modes)
xc_sampl = [-p/2;linspace(p/2,l-p/2,3)']; %sample locations of center thermal pulse

T_sampl = zeros(nNodes,length(xc_sampl)); % T nodal distribution samples
xc_sampl = [-p/2;0;0.01;0.05;0.085];

%generate corresponding temperature profiles
for ii = 1:length(xc_sampl)
    T_sampl(:,ii) = T_dyn_xc(xc_sampl(ii));
end

VMs_at_eq = 1; %do you want to compute VMs at thermal equilibrium? (set to 1 if yes, otherwise set to 0)
[VMs,static_eq,omega] = VM_thermal(BeamAssembly,T_sampl,VMs_at_eq,4);

%plot VMs__________________________________________________________________
T_sampl_2_plot = 1:length(xc_sampl);
mode2plot = [1,2,3,4];

%T distribution 
figure('units','normalized','position',[.3 .3 .4 .4]); hold on;
title('T distribution samples'); xlabel('x'); ylabel('T');
for jj = 1:length(T_sampl_2_plot)
    plot(nodes_x,T_sampl(:,T_sampl_2_plot(jj)))
end

%VMs
for ii = 1:length(mode2plot)
    
    VM2pl = mode2plot(ii);
    figure('units','normalized','position',[.3 .3 .4 .4]); hold on;
    title(['VM ',num2str(mode2plot(ii))]); xlabel('x');
    
    for jj = 1:length(T_sampl_2_plot)
        %decode output
        VM_iijj = decodeDofsNodes(VMs{jj}(:,VM2pl),nNodes,nDofPerNode);
        v_ii = VM_iijj(:,2); %tranverse displ.
        u_ii = VM_iijj(:,1); %long. displ.
        plot(nodes_x + u_ii, nodes_y + v_ii,'.-','markersize',10)
    end

end

%Equilibrium position
figure('units','normalized','position',[.3 .3 .4 .4]); hold on;
for ii = 1:length(T_sampl_2_plot)
   title('Static thermal equilibrium '); xlabel('x');
   ueq_ii = reshape(static_eq(:,ii),nDofPerNode,[]).';
   subplot 311; hold on;
   plot(nodes_x,ueq_ii(:,1)); 
   xlabel('x'); ylabel('u');
end
for ii = 1:length(T_sampl_2_plot)
   ueq_ii = reshape(static_eq(:,ii),nDofPerNode,[]).';
   subplot 312; hold on;
   plot(nodes_x,ueq_ii(:,2)); 
   xlabel('x'); ylabel('v');
end
scl = 50; %scaling
for ii = 1:length(T_sampl_2_plot)
   ueq_ii = reshape(static_eq(:,ii),nDofPerNode,[]).';
   subplot 313; hold on;
   plot(nodes_x + scl*ueq_ii(:,1), nodes_y + scl*ueq_ii(:,2),'*-'); 
   xlabel('x'); ylabel('deformation')
end


%% Reduced Order Model Construction 
% define number of equidistant samples
n_ROMs = 20;

xc_sampl = linspace(-p/2,l+p/2,n_ROMs)'; %sample locations of center thermal pulse

x0_sampl = xc_sampl - p/2; %left extreme of pulse
T_sampl = zeros(nNodes,n_ROMs); % T nodal distribution samples

%generate corresponding temperature profiles
for ii = 1:n_ROMs
    T_sampl(:,ii) = T_dyn_xc(xc_sampl(ii));
end

%construct model___________________________________________________________
% define how many VMs to use in modal reduction (also the corresponding
% modal derivatives are added)
number_VMs = 5;

ROMs.fullAssembly = BeamAssembly;
ROMs.parameters = xc_sampl'; 
ROMs.models = multiple_ROMs_thermal(BeamAssembly, T_sampl, number_VMs); 
% cell array, with in each cell a Model. Each model is a struct array with field V (for now, then could be extended to M, K, K2, K3)

%% Integrate nonlinear ROM

% settings for integration
tmax = T_th/4; % integration interval
h = T_fext/50; % time step for integration

% Initial condition: equilibrium
n_dof_r = size(ROMs.models{1}.V,2);

q0 = zeros(n_dof_r,1);
qd0 = zeros(n_dof_r,1);
qdd0 = zeros(n_dof_r,1);

% Instantiate object for nonlinear time integration
TI_NL_r = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear Residual evaluation function handle
residual_r = @(q,qd,qdd,t)residual_ROM_thermal(q,qd,qdd,t,ROMs,F_ext,xc_dyn,T_dyn_xc);

% integrate equations with Newmark scheme
TI_NL_r.Integrate(q0,qd0,qdd0,tmax,residual_r);

% postprocess output 
u_dyn_r = reduced_to_full_thermal(TI_NL_r.Solution.q,TI_NL_r.Solution.time,ROMs,xc_dyn,1:nDofsC);

u_dyn_r = BeamAssembly.unconstrain_vector(u_dyn_r); 
dyn_r.nlin.disp = decodeDofsNodes(u_dyn_r,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn_r.nlin.time = TI_NL_r.Solution.time;


%% Integrate nonlinear ROM (constant basis)

% The ROM with constant basis, can be obtained (neglecting computational
% efficiency) by setting all the basis of the different precomputed ROMs to
% the basis of the cold structure (the first one). In this way the rest of
% the code remains the same.
ROMs_cbasis = ROMs;
V_cold = ROMs.models{1}.V;
for ii = 1:n_ROMs
    ROMs_cbasis.models{ii}.V = V_cold;
end

% settings for integration
tmax = T_th/4; % integration interval
h = T_fext/50; % time step for integration

% Initial condition: equilibrium
n_dof_r = size(ROMs.models{1}.V,2);

q0 = zeros(n_dof_r,1);
qd0 = zeros(n_dof_r,1);
qdd0 = zeros(n_dof_r,1);

% Instantiate object for nonlinear time integration
TI_NL_r_cbasis = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear Residual evaluation function handle
residual_r_cbasis = @(q,qd,qdd,t)residual_ROM_thermal(q,qd,qdd,t,ROMs_cbasis,F_ext,xc_dyn,T_dyn_xc);

% integrate equations with Newmark scheme
TI_NL_r_cbasis.Integrate(q0,qd0,qdd0,tmax,residual_r_cbasis);

% postprocess output 
u_dyn_r_cbasis = reduced_to_full_thermal(TI_NL_r_cbasis.Solution.q,TI_NL_r_cbasis.Solution.time,ROMs_cbasis,xc_dyn,1:nDofsC);

u_dyn_r_cbasis = BeamAssembly.unconstrain_vector(u_dyn_r_cbasis); 
dyn_r_cb.nlin.disp = decodeDofsNodes(u_dyn_r_cbasis,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn_r_cb.nlin.time = TI_NL_r_cbasis.Solution.time;

%% Plot comparison
% plot results of nonlinear dynamic analysis_______________________________
plot_dof_location = 0.3; %percentage of beam length
node2plot = find_node(plot_dof_location*l,0,[],Nodes); % node where to put the force
plot_dof = get_index(node2plot, BeamAssembly.Mesh.nDOFPerNode );

figure('units','normalized','position',[.1 .1 .8 .8]); 
subplot 311; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,1,:)),'-','color','r');
plot(dyn_r.nlin.time,squeeze(dyn_r.nlin.disp(node2plot,1,:)),'--','color','b');
plot(dyn_r_cb.nlin.time,squeeze(dyn_r_cb.nlin.disp(node2plot,1,:)),'-.','color','m');
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,1,:)),'--','color','cyan');
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,1,:)),'--','color','g');
xlabel('t [s]'); ylabel('u [m]'); grid on; axis tight; legend('full','ROM','ROM const V','lin','quasistatic');
title('comparison of ROM and full solutions')

subplot 312; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,2,:)),'-','color','r');
plot(dyn_r.nlin.time,squeeze(dyn_r.nlin.disp(node2plot,2,:)),'--','color','b');
plot(dyn_r_cb.nlin.time,squeeze(dyn_r_cb.nlin.disp(node2plot,2,:)),'-.','color','m');
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,2,:)),'--','color','cyan');
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,2,:)),'--','color','g');
xlabel('t [s]'); ylabel('v [m]'); grid on; axis tight; legend('full','ROM','ROM const V','lin','quasistatic');

subplot 313; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,3,:)),'-','color','r');
plot(dyn_r.nlin.time,squeeze(dyn_r.nlin.disp(node2plot,3,:)),'--','color','b');
plot(dyn_r_cb.nlin.time,squeeze(dyn_r_cb.nlin.disp(node2plot,3,:)),'-.','color','m');
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,3,:)),'--','color','cyan');
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,3,:)),'--','color','g');
xlabel('t [s]'); ylabel('w [-]'); grid on; axis tight; legend('full','ROM','ROM const V','lin','quasistatic');

return
%% plot comparison

% for ii = 1:size(TI_NL_r.Solution.q,1)
%     figure
%     plot(TI_NL_r.Solution.time,TI_NL_r.Solution.q(ii,:))
% end

% postprocess output 
u_dyn_r = reduced_to_full_thermal(TI_NL_r.Solution.q,TI_NL_r.Solution.time,ROMs,xc_dyn,1:nDofsC);

u_dyn_r = BeamAssembly.unconstrain_vector(u_dyn_r); 
dyn_r.nlin.disp = decodeDofsNodes(u_dyn_r,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn_r.nlin.time = TI_NL_r.Solution.time;

%quasistatic response
u_qs = quasistatic_from_ROMs(TI_NL_r.Solution.time,ROMs,xc_dyn);
u_qs = BeamAssembly.unconstrain_vector(u_qs);
u_qs = decodeDofsNodes(u_qs,nNodes,nDofPerNode);

% plot results of nonlinear dynamic analysis_______________________________

plot_dof_location = 0.7; %percentage of length of beam
node2plot = find_node(plot_dof_location*l,0,[],Nodes); % node where to put the force
plot_dof = get_index(node2plot, BeamAssembly.Mesh.nDOFPerNode );

figure;
subplot 311; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,1,:)),'-');
plot(dyn_r.nlin.time,squeeze(dyn_r.nlin.disp(node2plot,1,:)),'x-');
plot(dyn_r.nlin.time,squeeze(u_qs(node2plot,1,:)),'o-');
xlabel('t'); ylabel('u'); grid on; axis tight; 

subplot 312; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,2,:)),'-');
plot(dyn_r.nlin.time,squeeze(dyn_r.nlin.disp(node2plot,2,:)),'x-');
plot(dyn_r.nlin.time,squeeze(u_qs(node2plot,2,:)),'o-');
xlabel('t'); ylabel('v'); grid on; axis tight; 

subplot 313; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,3,:)),'-');
plot(dyn_r.nlin.time,squeeze(dyn_r.nlin.disp(node2plot,3,:)),'x-');
plot(dyn_r.nlin.time,squeeze(u_qs(node2plot,3,:)),'o-');
xlabel('t'); ylabel('w'); grid on; axis tight; 






