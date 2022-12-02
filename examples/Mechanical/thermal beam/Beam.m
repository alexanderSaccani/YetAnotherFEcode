% ALEXANDER SACCANI - PHD CANDIDATE, ETH ZURICH, 
% created: 1.12.2022
% last modified: 1.12.2022

% Straight temperature dependent beam template
% single campled, axial loading

clearvars ; close all; clc

%% Model construction

%GEOMETRY__________________________________________________________________
l = 0.2; 
h = 1e-3;
b = 1e-2; 


%MATERIAL__________________________________________________________________
E       = 200e9;  % 70e9 % 200e9 % Young's modulus
rho     = 7850; % 2700 % 7850 % density
nu      = 0.3;    % nu
kappa   = 1e8; % material damping modulus 1e8
alpha_T = 11.7e-6; % thermal expansion coefficient

myThermalBeamMaterial = KirchoffMaterial();
set(myThermalBeamMaterial,'YOUNGS_MODULUS', E,'DENSITY',rho,...
    'POISSONS_RATIO',nu,'THERMAL_EXPANSION_COEFFICIENT', alpha_T);


%MESH______________________________________________________________________
nElements = 10;

%define nodes
dx = l/nElements; %x spacing between nodes

nodes_x = (0:dx:l).'; %x coordinates for nodes
nodes_y = zeros(length(nodes_x),1); %y coordinates for nodes

nNodes = length(nodes_x); %number of nodes
nDofPerNode = 3;

%create mesh from nodes
Nodes = [nodes_x, nodes_y];
BeamMesh = Mesh(Nodes); 

%define elements 
Elements = [1:nNodes-1;2:nNodes].'; % ROW -> # element, 
                                    % COLUMN -> # of node connected
myElementConstructor = @()ThermalBeamElement(b, h, myThermalBeamMaterial); %type of element

%create element table from elements
BeamMesh.create_elements_table(Elements,myElementConstructor);


%BOUNDARY CONDITIONS_______________________________________________________

%BeamMesh.set_essential_boundary_condition([1 BeamMesh.nNodes],[1 2 3],0) % Doubly clamped beam
BeamMesh.set_essential_boundary_condition(1,[1 2 3],0) % clamped beam  (node,[constrained dofs],imposed disp.)

BeamAssembly = Assembly(BeamMesh);

nDofsF = BeamAssembly.Mesh.EBC.nDOFs; % #dofs free (without bc)
nDofsC = nDofsF - size(BeamAssembly.Mesh.EBC.constrainedDOFs,1); % #dofs constrained (with bc)


%GET USEFUL QUANTITIES FROM THE CREATED MODEL______________________________
K = BeamAssembly.stiffness_matrix(); % should be tangent stiffness matrix for T = 0, u = 0
M = BeamAssembly.mass_matrix();
C = BeamAssembly.damping_matrix();

% u0 = zeros(nDofsF,1);
% T0 = 0*ones(nNodes,1);
% [K0,~] = ThermalBeamAssembly.tangent_stiffness_and_force(u0,T0);
% 
% diff = (K - K0);
% diff = full(diff);

K_C = BeamAssembly.constrain_matrix(K); % should be tangent stiffness matrix for T = 0, u = 0
M_C = BeamAssembly.constrain_matrix(M);

%% Modal Analysis


n_VMs = 3; % first n_VMs modes with lowest frequency calculated 

[VM,omega2] = eigs(K_C,M_C,n_VMs,'SM'); 
omega = sqrt(diag(omega2));
%normalize with respect to max. amplitude
for ii = 1:n_VMs
    VM(:,ii) = VM(:,ii)/max(sqrt(sum(VM(:,ii).^2,2)));
end
VM = BeamAssembly.unconstrain_vector(VM); %vibration modes

%decode output
VM = decodeDofsNodes(VM,nNodes,nDofPerNode);

%plot VMs
mode2plot = [1,2,3];
for ii = 1:length(mode2plot)
  
    VM_ii = squeeze(VM(:,:,mode2plot(ii)));
    
    u_ii = VM_ii(:,1); % long. displ.
    v_ii = VM_ii(:,2); % transv. displ
    r_ii = VM_ii(:,3); % rotation
    
    figure
    subplot 311
    plot(nodes_x,nodes_y,'.-','markersize',10); hold on;
    plot(nodes_x + u_ii, nodes_y + v_ii,'.-','markersize',10)
    title(['VM ',num2str(mode2plot(ii)),', Frequency = '...
    num2str(omega(mode2plot(ii))/(2*pi)) ' Hz'])
    
    subplot 312
    plot(nodes_x,v_ii,'-'); title('tranverse disp.'); xlabel('x')
    subplot 313
    plot(nodes_x,u_ii,'-'); title('longitudinal disp.'); xlabel('x')
    
    
end


%% define Thermal modes

% thermal mode 1 (uniform temp. distribution)
th1 = ones(nNodes,1); 

% thermal mode 2 (uniform temp. distribution over one half of the beam)
th2 = zeros(nNodes,1);
th2(1:round(nNodes/2)) = ones(round(nNodes/2),1); 


%% Static analysis 

% Nodal force
F = zeros(nDofsF,1);

% transverse load at midspan
nodeF = find_node(l/2,0,[],Nodes); % node where to put the force
node_force_dofs = get_index(nodeF, nDofPerNode );
F(node_force_dofs(2)) = 0; 

% Temperature field for static analysis
T_st = 100*th1; % static temperature field (half beam is hot, half cold)

% % linear static solution (wrong since the temperature is not
% % considered!)          PLEASE CHECK THIS!!!!
%  u_lin = ThermalBeamAssembly.solve_system(K, F);     %this is wrong since the temperature must enter here!
% static.lin = reshape(u_lin,3,[]).';	% Linear response

% nonlinear static solution (here instead the temperature is considered)
[K,~] = BeamAssembly.tangent_stiffness_and_force(zeros(nDofsF,1),zeros(nNodes));
BeamAssembly.DATA.K = K;    % this value of K is used only for first guess solution (T not considered in the linear guess!)
[~,u] = static_equilibrium_thermal(BeamAssembly, F, T_st,'display', 'iter-detailed');
static.nlin = reshape(u,nDofPerNode,[]).';     % Nonlinear response

%plot static displacements
ud = static.nlin(:,1);
vd = static.nlin(:,2);
wd = static.nlin(:,3);
scl = 50;   %scaling factor

figure
subplot 311; plot(nodes_x,ud); xlabel('x'); ylabel('u'); title('long. displ');
subplot 312; plot(nodes_x,vd); xlabel('x'); ylabel('v'); title('tranv. displ');
subplot 313; plot(nodes_x,nodes_y,'.-','markersize',10,'color','b'); hold on;
plot(nodes_x + scl*ud, nodes_y + scl*vd, '.-', 'markersize',10,'color','r')
title('static deformation shape'), xlabel('x')


%% Dynamic response using Implicit Newmark

%force and thermal load____________________________________________________

%define external forcing
om_fext = (omega(1) + omega(2))/2;
T_fext = 2*pi/om_fext; % period of ex. forcing
F_ampl = 100; 

node_fext = find_node(l/2,0,[],Nodes); % node where to put the force
fext_dof = get_index(node_fext, nDofPerNode );

%construct forcing vector
F = zeros(nDofsF,1);
F(fext_dof(1)) = F_ampl;
F_ext = @(t) F*sin(om_fext*t);

%define thermal load 
eps = 0.05; % om_fex/om_th (th = thermal forcing)
om_th = eps*om_fext;  % angular frequency of thermal mode oscillation in time
T_th = 2*pi/om_th; % period of temp. oscillations

T_dyn = @(t) 100*th1*sin(om_th*t); % amplitude of thermal mode

% Integrate nonlinear EOM__________________________________________________

% settings for integration
tmax = 2*T_th; % integration interval
h = T_fext/50; % time step for integration

% Initial condition: equilibrium
q0 = zeros(nDofsC,1);
qd0 = zeros(nDofsC,1);
qdd0 = zeros(nDofsC,1);

% Precompute data for Assembly object
BeamAssembly.DATA.M = M;
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.C = C;

% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_thermal_nonlinear(q,qd,qdd,t,...
    BeamAssembly,F_ext, T_dyn);

% integrate equations with Newmark scheme
TI_NL.Integrate(q0,qd0,qdd0,tmax,residual);

% store solution in more suitable array (only displacements)
u_dyn = BeamAssembly.unconstrain_vector(TI_NL.Solution.q); 
dyn.nlin.disp = decodeDofsNodes(u_dyn,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn.nlin.time = TI_NL.Solution.time;

% plot results of nonlinear dynamic analysis_______________________________

plot_dof_location = 0.3; %percentage of length of beam
node2plot = find_node(plot_dof_location*l,0,[],Nodes); % node where to put the force
plot_dof = get_index(node2plot, BeamAssembly.Mesh.nDOFPerNode );

figure;
subplot 311
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,1,:)),'-');
xlabel('t'); ylabel('u'); grid on; axis tight; 

subplot 312
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,2,:)),'-');
xlabel('t'); ylabel('v'); grid on; axis tight; 

subplot 313
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,3,:)),'-');
xlabel('t'); ylabel('w'); grid on; axis tight; 










