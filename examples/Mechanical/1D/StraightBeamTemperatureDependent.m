%% A temperature-dependent beam model

clear all; close all; clc
%% parameters
% geometry
l = 0.2; %1 % 0.2
h = 1e-3;
b = 1e-2; %1e-1 % 1e-2

w   = 5; %1/b*l is the width of the sin^2 graph of temperature pulse on the beam
c = 0.215; % location of the pulse center: factor of beam length
n_VMs = 5; % number of vibration modes
n_BMs = 5; % number of buckling modes
x_c = 0.08:0.04:0.92; % Center of the sin^2: factor of beam length

% mesh
nElements = 3;
dx = l/nElements;
outdof = 3*floor(nElements/4) + 2; % output degree of freedom for displaying results

% material properties
E       = 200e9;  % 70e9 % 200e9 % Young's modulus
rho     = 7850; % 2700 % 7850 % density
nu      = 0.3;    % nu
kappa   = 1e8; % material damping modulus 1e8
alpha_T = 11.7e-6; % thermal expansion coefficient
loadfactor = 1; % mechanical load factor 
%% Structural
% Material
myBeamMaterial  = KirchoffMaterial();
set(myBeamMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'DAMPING_MODULUS',kappa);
D               = myBeamMaterial.get_stress_strain_matrix_2D();

% Element
myElementConstructor = @()BeamElement(b, h, myBeamMaterial); % same element all across the domain

% Mesh
x = (0:dx:l).';
nNodes = size(x,1);
Nodes = [x, zeros(nNodes,1)];

Elements = [1:nNodes-1;2:nNodes].';
BeamMesh = Mesh(Nodes);
BeamMesh.create_elements_table(Elements,myElementConstructor);

% Assembly
u0 = zeros(nNodes*BeamMesh.nDOFPerNode,1);
BeamAssembly = Assembly(BeamMesh);
[K, F] = BeamAssembly.tangent_stiffness_and_force(u0);
M = BeamAssembly.mass_matrix();

u0 = randi(5,BeamAssembly.Mesh.nDOFs,1);
F2 = BeamAssembly.vector('F2',u0,u0);
T2 = BeamAssembly.tensor('T2',[BeamMesh.nDOFs, BeamMesh.nDOFs, BeamMesh.nDOFs], [2,3]);
F2check = ttv(T2,{u0,u0},[2,3]);
norm(F2check.data - F2)/norm(F2)

F3 = BeamAssembly.vector('F3',u0,u0,u0);
T3 = BeamAssembly.tensor('T3',[BeamMesh.nDOFs, BeamMesh.nDOFs, BeamMesh.nDOFs, BeamMesh.nDOFs], [2,3,4]);
F3check = ttv(T3,{u0,u0,u0},[2,3,4]);
norm(F3check.data - F3)/norm(F3)


[~,F] = BeamAssembly.tangent_stiffness_and_force(u0);
norm(K*u0 + F2 + F3 - F)/norm(F)

S = BeamAssembly.scalar('strain_energy',u0);

%% Thermal beam
% Material
myThermalBeamMaterial = KirchoffMaterial();
set(myThermalBeamMaterial,'YOUNGS_MODULUS', E,'DENSITY',rho,'POISSONS_RATIO',nu,'THERMAL_EXPANSION_COEFFICIENT', alpha_T);

% Element
myElementConstructor = @()ThermalBeamElement(b, h, myThermalBeamMaterial); % same element all across the domain

% Mesh
ThermalBeamMesh = Mesh(Nodes);
ThermalBeamMesh.create_elements_table(Elements,myElementConstructor);
nNodes = ThermalBeamMesh.nNodes;

nodes_x = ThermalBeamMesh.nodes(:,1);
nodes_y = ThermalBeamMesh.nodes(:,2);

% Set dirichlet DOFs
%ThermalBeamMesh.set_essential_boundary_condition([1 ThermalBeamMesh.nNodes],[1 2 3],0) % Doubly clamped beam
ThermalBeamMesh.set_essential_boundary_condition([1],[1 2 3],0) % clamped beam
ThermalBeamAssembly = Assembly(ThermalBeamMesh);

nDofsF = ThermalBeamAssembly.Mesh.EBC.nDOFs;
nDofsC = nDofsF - size(ThermalBeamAssembly.Mesh.EBC.constrainedDOFs,1);
nNodes = ThermalBeamAssembly.Mesh.nNodes;

% force
% uniform transverse external force
%f = ThermalBeamAssembly.uniform_body_force();

%% eigenvalue analysis

K = ThermalBeamAssembly.stiffness_matrix(); % should be the stiffness matrix for T = 0
M = ThermalBeamAssembly.mass_matrix();
C = ThermalBeamAssembly.damping_matrix();

% u0 = zeros(nDofsF,1);
% T0 = 0*ones(nNodes,1);
% [K0,~] = ThermalBeamAssembly.tangent_stiffness_and_force(u0,T0);
% 
% diff = (K - K0);
% diff = full(diff);

n_VMs = 2; % first n_VMs modes with lowest frequency calculated 
K_red = ThermalBeamAssembly.constrain_matrix(K); % at which temperature
M_red = ThermalBeamAssembly.constrain_matrix(M);
[VM,omega2] = eigs(K_red,M_red,n_VMs,'SM'); 
omega = sqrt(diag(omega2));
%normalize with respect to max. amplitude
for ii = 1:n_VMs
    VM(:,ii) = VM(:,ii)/max(sqrt(sum(VM(:,ii).^2,2)));
end
VM = ThermalBeamAssembly.unconstrain_vector(VM);

% %plot VMs
% mode2plot = [1,2,3];
% for ii = 1:length(mode2plot)
%     v2plot = reshape(VM(:,mode2plot(ii)),3,[]);
%     %plot VMs
%     ud = v2plot(1,:)';
%     vd = v2plot(2,:)';
%     rot = v2plot(3,:)'; %this is the rotation (here we do not plot it)
%     figure
%     subplot 311
%     title(['VM ',num2str(mode2plot(ii)),' Frequency = ' num2str(omega(mode2plot(ii))/(2*pi)) ' Hz'])
%     plot(nodes_x,nodes_y,'-.','markersize',10); hold on;
%     plot(nodes_x + ud, nodes_y + vd)
%     subplot 312
%     plot(nodes_x,vd,'-'); title('tranverse disp.');
%     subplot 313
%     plot(nodes_x,ud,'-'); title('longitudinal disp.');
% end

%% define Thermal modes

% thermal mode 1 (uniform temp. distribution)
th1 = ones(nNodes,1); 

% thermal mode 2 (uniform temp. distribution over one half of the beam)
th2 = zeros(nNodes,1);
th2(1:round(nNodes/2)) = ones(round(nNodes/2),1); %


%% Static analysis 

% Nodal force
F = zeros(nDofsF,1);
node2plot = find_node(l/2,0,[],Nodes); % node where to put the force
node_force_dofs = get_index(node2plot, ThermalBeamAssembly.Mesh.nDOFPerNode );
F(node_force_dofs(2)) = 0;

% Temperature field for static analysis
T_st = 100*th2; % static temperature field (half beam is hot, half cold)

% % linear static solution (wrong since the temperature is not
% % considered!)          PLEASE CHECK THIS!!!!
%  u_lin = ThermalBeamAssembly.solve_system(K, F);                     %this is wrong since the temperature must enter here!
% static.lin = reshape(u_lin,3,[]).';	% Linear response

% nonlinear static solution (here instead the temperature is considered)
[K,~] = ThermalBeamAssembly.tangent_stiffness_and_force(zeros(nDofsF,1),zeros(nNodes));
ThermalBeamAssembly.DATA.K = K;    % this value of K is used only for first guess solution (T not considered in the linear guess!)
[~,u] = static_equilibrium_thermal(ThermalBeamAssembly, F, T_st,'display', 'iter-detailed');
static.nlin = reshape(u,3,[]).';     % Nonlinear response

%plot static displacements
ud = static.nlin(:,1);
vd = static.nlin(:,2);
wd = static.nlin(:,3);
scl = 50;   %scaling factor

% figure
% subplot 311; plot(nodes_x,ud);
% subplot 312; plot(nodes_x,vd);
% subplot 313; plot(nodes_x,nodes_y,'.-','markersize',10,'color','b'); hold on;
% plot(nodes_x + scl*ud, nodes_y + scl*vd, '.-', 'markersize',10,'color','r')
% title('static deformation'), xlabel('x')


%% Dynamic response using Implicit Newmark

%thermal load
eps = 0.05;%1e-2; 
om_th = eps*omega(1);  % angular frequency of thermal mode oscillation in time
T_temp = 2*pi/om_th; % period of temp. oscillations
T_dyn = @(t) 100*th1*sin(om_th*t); % amplitude of thermal mode

%external forcing
om_forcing = (omega(1) + omega(2))/2;
F_ampl = 100; 

F = zeros(nDofsF,1);
node2plot = find_node(l/2,0,[],Nodes); % node where to put the force
node_force_dofs = get_index(node2plot, ThermalBeamAssembly.Mesh.nDOFPerNode );
F(node_force_dofs(1)) = F_ampl;

F_ext = @(t) F*sin(om_forcing*t);

T_forc =  2*pi/om_forcing; % time period of forcing

% Initial condition: equilibrium
q0 = zeros(nDofsC,1);
qd0 = zeros(nDofsC,1);
qdd0 = zeros(nDofsC,1);

% time step for integration
h = T_forc/50;

% Precompute data for Assembly object
ThermalBeamAssembly.DATA.M = M;
ThermalBeamAssembly.DATA.K = K;
ThermalBeamAssembly.DATA.C = C;

%% Integration of linear EOMs. Does it have sens?
% % Instantiate object for linear time integration
% TI_lin = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
% 
% % Linear Residual evaluation function handle
% residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,PlateAssembly,F_ext);
% 
% % Linearized Time Integration
% tmax = 10*T; 
% TI_lin.Integrate(q0,qd0,qdd0,tmax,residual_lin);
% 
% % obtain full solution
% TI_lin.Solution.u = PlateAssembly.unconstrain_vector(TI_lin.Solution.q);
% 
% % Animate solution on Mesh (very slow)
% % AnimateFieldonDeformedMesh(myMesh.Nodes,myMesh.Elements,TI_lin.Solution.u ,'factor',1,'index',1:3,'filename','lineardisp')

%% Integration of nonlinear EOMs
% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

% Linear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_thermal_nonlinear(q,qd,qdd,t,...
    ThermalBeamAssembly,F_ext, T_dyn);

% Nonlinear Time Integration
tmax = 50*T_forc; 
TI_NL.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL.Solution.u = ThermalBeamAssembly.unconstrain_vector(TI_NL.Solution.q);

%% save displacements in another format

n_tsamp = length(TI_NL.Solution.time);
nDOFPerNode = ThermalBeamAssembly.Mesh.nDOFPerNode;

dyn.nlin.disp = zeros(nNodes,nDOFPerNode,n_tsamp); % (node, dof of node, tsamp)
u = ThermalBeamAssembly.unconstrain_vector(TI_NL.Solution.q); %unconstrained vector of displacements
for ii = 1:n_tsamp
   disp_ii = u(:,ii);
   dyn.nlin.disp(:,:,ii) = (reshape(disp_ii,nDOFPerNode,[])).';
end
dyn.nlin.time = TI_NL.Solution.time;


%% plot displacement of midspan
plot_dof_location = 0.3; %percentage of length of beam
node2plot = find_node(plot_dof_location*l,0,[],Nodes); % node where to put the force
plot_dof = get_index(node2plot, ThermalBeamAssembly.Mesh.nDOFPerNode );

figure;
subplot 311
%plot(TI_NL.Solution.time, TI_NL.Solution.u(plot_dof(1),:),'DisplayName','v midspan'); hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,1,:)),'-');
xlabel('t'); ylabel('u'); grid on; axis tight; 

subplot 312
%plot(TI_NL.Solution.time, TI_NL.Solution.u(plot_dof(2),:),'DisplayName','u midspan'); hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,2,:)),'-');
xlabel('t'); ylabel('v'); grid on; axis tight; 

subplot 313
%plot(TI_NL.Solution.time, TI_NL.Solution.u(plot_dof(3),:),'DisplayName','w midspan'); hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,3,:)),'-');
xlabel('t'); ylabel('w'); grid on; axis tight; 










