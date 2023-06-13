close all
clc
clearvars

%% Model construction

%MAIN BEAM

%GEOMETRY__________________________________________________________________
l = 0.1; 
h = 1e-3;
b = 1e-2; 


%MATERIAL__________________________________________________________________
E       = 200e9;  % 70e9 % 200e9 % Young's modulus
rho     = 7850; % 2700 % 7850 % density
nu      = 0.3;    % nu
kappa   = 1e8; % material damping modulus 1e8
alpha_T = 11.7e-6; % thermal expansion coefficient

BeamMaterial = KirchoffMaterial();
set(BeamMaterial,'YOUNGS_MODULUS', E,'DENSITY',rho,...
    'POISSONS_RATIO',nu,'THERMAL_EXPANSION_COEFFICIENT', alpha_T);


%MESH______________________________________________________________________
nElements = 20;

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
myElementConstructor = @()BeamElement(b, h, BeamMaterial); %type of element

%create element table from elements
BeamMesh.create_elements_table(Elements,myElementConstructor);


%BOUNDARY CONDITIONS_______________________________________________________

%BeamMesh.set_essential_boundary_condition([1 BeamMesh.nNodes],[1 2 3],0) % Doubly clamped beam
BeamMesh.set_essential_boundary_condition(1,[1 2 3],0) % clamped beam  (node,[constrained dofs],imposed disp.)
BeamMesh.set_essential_boundary_condition(size(Nodes,1),[1],0)

BeamAssembly = Assembly(BeamMesh);

nDofsF = BeamAssembly.Mesh.EBC.nDOFs; % #dofs free (without bc)
nDofsC = nDofsF - size(BeamAssembly.Mesh.EBC.constrainedDOFs,1); % #dofs constrained (with bc)


% ADDED SPRING
I = 1/12*h^3*b;
k_beam = 3*E*I/l^3;

% r = 1;
% 
% k_spring = r*k_beam;


node_spring = find_node(l,0,[],Nodes); % node where to put the force
spring_dof = get_index(node_spring, nDofPerNode );

n_dof_spring = spring_dof(2);

% n_dof_spring = (size(Nodes,1)-1)*3+2;

%GET USEFUL QUANTITIES FROM THE CREATED MODEL______________________________
M = BeamAssembly.mass_matrix();
C = BeamAssembly.damping_matrix();


M_C = BeamAssembly.constrain_matrix(M);


%% Modal Analysis varying spring stiffness

r_all = [0,2,5,10,20,50].';

n_VMs = 2; % first n_VMs modes with lowest frequency calculated 

VM = cell(length(r_all));
omega = zeros(n_VMs,length(r_all));

for ii = 1:length(r_all)
    
    k_spring = r_all(ii)*k_beam;
    [K,~] = Ktg_force(BeamAssembly,k_spring,n_dof_spring,zeros(3*size(Nodes,1),1)); % tangent stiffness matrix for u = 0
    
    K_C_ii = BeamAssembly.constrain_matrix(K); %tangent stiffness matrix for u = 0
    
    [VM_ii,omega2_ii] = eigs(K_C_ii,M_C,n_VMs,'SM'); 
    omega_ii = sqrt(diag(omega2_ii));
    %normalize with respect to max. amplitude
    for jj = 1:n_VMs
        VM_ii(:,jj) = VM_ii(:,jj)/vecnorm(VM_ii(:,jj));
    end
    VM_ii = BeamAssembly.unconstrain_vector(VM_ii); %vibration modes
    
    VM{ii} = VM_ii;
    omega(:,ii) = omega_ii.';  
    
end

color_list = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE',...
    '#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};

%plot VMs
scale_factor = 1;
mode2plot = 1:n_VMs;
for ii = 1:length(mode2plot)
    
    VM2pl = mode2plot(ii);
    figure('units','normalized','position',[.3 .3 .4 .4]); hold on;
    title(['VM ',num2str(mode2plot(ii))]);
    
    for jj = 1:length(r_all)
        %decode output
        VM_iijj = decodeDofsNodes(VM{jj}(:,VM2pl),nNodes,nDofPerNode);
        PlotFieldonDeformedMesh(Nodes,Elements,[VM_iijj(:,1),VM_iijj(:,2)],'factor',scale_factor, 'color',color_list{jj});
        PlotFieldonDeformedMesh(Nodes,Elements,[-VM_iijj(:,1),-VM_iijj(:,2)],'factor',scale_factor, 'color',color_list{jj});
    end

end

disp('natural frequencies are')
disp(omega);


%% define spring stiffness variation in time 

%parametric variation of k(r)
k_r = @(r) k_beam*r; 

%temporal variation of r(t)
eps = 0.01; 
om_k = eps*min(min(omega));
T_k = 2*pi/om_k;

% %sine variation
% A_max = 20;
% r_0 = 0;
% r_t = @(t) r_0 + A_max*sin(om_k*t);

%ramp variation
A_max = 50;
r_t = @(t) A_max/T_k*4*t;


k_t = @(t) k_r(r_t(t));


%% define forcing

%define external forcing
om_fext = omega(1,1) + 0.35*(omega(2,1) - omega(1,1)); % (omega(1,1) + omega(2,1))/2;
T_fext = 2*pi/om_fext ; % period of ex. forcing
F_ampl = 10; 

node_fext = find_node(l/2,0,[],Nodes); % node where to put the force
fext_dof = get_index(node_fext, nDofPerNode );

%construct forcing vector
F = zeros(nDofsF,1);
F(fext_dof(2)) = F_ampl;
F_ext = @(t) F*sin(om_fext*t);

%% nonlinear time integration

% Integrate nonlinear EOM__________________________________________________

% settings for integration
tint = T_k/4; % integration interval
h = T_fext/30; % time step for integration

% Initial condition: equilibrium
q0 = zeros(nDofsC,1);
qd0 = zeros(nDofsC,1);
qdd0 = zeros(nDofsC,1);

% Precompute data for Assembly object
BeamAssembly.DATA.M = M;
[K_0,~] = Ktg_force(BeamAssembly,0,n_dof_spring,zeros(3*size(Nodes,1),1));
BeamAssembly.DATA.K = K_0;
BeamAssembly.DATA.C = C;

% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear Residual evaluation function handle
residual = @(q,qd,qdd,t) residual_beam_spring(q,qd,qdd,t,BeamAssembly,...
                        F_ext,k_t,n_dof_spring);

% integrate equations with Newmark scheme
TI_NL.Integrate(q0,qd0,qdd0,tint,residual);

% store solution in more suitable array (only displacements)
u_dyn = BeamAssembly.unconstrain_vector(TI_NL.Solution.q); 
ud_dyn = BeamAssembly.unconstrain_vector(TI_NL.Solution.qd); 
dyn.nlin.disp = decodeDofsNodes(u_dyn,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn.nlin.time = TI_NL.Solution.time;

dyn.nlin.vel = decodeDofsNodes(ud_dyn,nNodes,nDofPerNode);

% plot results of nonlinear dynamic analysis_______________________________

plot_dof_location = 2/3; %percentage of length of beam
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

figure
t = linspace(0,tint,1000);
plot(t,k_t(t));


%% reduced order model (constant basis)
r_samples = 0;
k_samples = r_samples*k_beam;
nVMs_ROM = 2;
logic_MD = 1;
[ROMS,k_ROMs] = construct_ROM(k_samples,BeamAssembly,nDofsF,n_dof_spring,nVMs_ROM,logic_MD);

V = ROMS{1}.V;

nDofs_red = size(V,2);

% Integrate nonlinear ROM__________________________________________________

% settings for integration
h = T_fext/50; % time step for integration

% Initial condition: equilibrium
q0 = zeros(nDofs_red,1);
qd0 = zeros(nDofs_red,1);
qdd0 = zeros(nDofs_red,1);


% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear red Residual evaluation function handle;                  
residual = @(q,qd,qdd,t) residual_beam_spring_ROM(q,qd,qdd,t,BeamAssembly,F_ext,k_t,n_dof_spring,V);

% integrate equations with Newmark scheme
TI_NL.Integrate(q0,qd0,qdd0,tint,residual);

% store solution in more suitable array (only displacements)
u_red = TI_NL.Solution.q;
ud_red = TI_NL.Solution.qd;

u_full = get_full_solution(u_red,V);
ud_full = get_full_solution(ud_red,V);

u_full = BeamAssembly.unconstrain_vector(u_full); 
ud_full = BeamAssembly.unconstrain_vector(ud_full); 

dyn.nlin_red.disp = decodeDofsNodes(u_full,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn.nlin_red.vel = decodeDofsNodes(ud_full,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn.nlin_red.time = TI_NL.Solution.time;


%%
% plot results of nonlinear dynamic analysis_______________________________

plot_dof_location = 0.4;%0.9%0.75;%2/3; %percentage of length of beam
node2plot = find_node(plot_dof_location*l,0,[],Nodes); % node where to put the force
plot_dof = get_index(node2plot, BeamAssembly.Mesh.nDOFPerNode );

figure;
subplot 311; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,1,:)),'-');
plot(dyn.nlin_red.time,squeeze(dyn.nlin_red.disp(node2plot,1,:)),'-');
xlabel('t'); ylabel('u'); grid on; axis tight; 

subplot 312; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,2,:)),'-');
plot(dyn.nlin_red.time,squeeze(dyn.nlin_red.disp(node2plot,2,:)),'-');
xlabel('t'); ylabel('v'); grid on; axis tight; 

subplot 313; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,3,:)),'-');
plot(dyn.nlin_red.time,squeeze(dyn.nlin_red.disp(node2plot,3,:)),'-');
xlabel('t'); ylabel('w'); grid on; axis tight; 

figure
t = linspace(0,tint,1000);
plot(t,k_t(t));


%% reduced order model (changing basis)

%construct local ROMs
r_samples = linspace(0,A_max,3);
k_samples = r_samples*k_beam;
nVMs_ROM = 2;
logic_MD = 1;
[ROMS,k_ROMs,k_validity_range] = construct_ROM(k_samples,BeamAssembly,nDofsF,n_dof_spring,nVMs_ROM,logic_MD);

%check subspace between basis
anglSub = 180/pi*subspace(ROMS{1}.V,ROMS{length(r_samples)}.V);

% Initial condition: equilibrium
u_f_F = zeros(nDofsC,1);
v_f_F = zeros(nDofsC,1);
a_f_F = zeros(nDofsC,1);

% tmin = 0;
% tmax = tint;
% u_full = [];
% time_red = [];
% u_full_red = [];
% 
% ii = 1;
% 
% while tmin < tmax
% 
% q0_ii = ROMS{ii}.V'*u_f_F;   
% qd0_ii = ROMS{ii}.V'*v_f_F; 
% qdd0_ii = ROMS{ii}.V'*a_f_F; 
%     
% term_cond = @(q,qd,t) exit_integration(q,qd,t,k_t,k_samples(ii+1));
% 
% % Instantiate object for nonlinear time integration
% TI_NL = ImplicitNewmark_mod('timestep',h,'alpha',0.005);
% 
% % nonlinear red Residual evaluation function handle;                  
% residual = @(q,qd,qdd,t) residual_beam_spring_ROM(q,qd,qdd,t,BeamAssembly,F_ext,k_t,n_dof_spring,ROMS{ii}.V);
% 
% % integrate equations with Newmark scheme
% TI_NL.Integrate(q0_ii,qd0_ii,qdd0_ii,tmin,tmax,residual,term_cond);
% 
% % store solution in more suitable array (only displacements)
% u_red = TI_NL.Solution.q;
% u_full_ii = get_full_solution(u_red,ROMS{ii}.V);
% 
% %save final displacements and velocities
% u_f_F = u_full_ii(:,end);
% v_f_F = ROMS{ii}.V*(TI_NL.Solution.qd(:,end));
% %accelerations ?? 
% 
% %save full disp and time vector
% u_full_ii = BeamAssembly.unconstrain_vector(u_full_ii); 
% u_full_red = [u_full_red,u_full_ii];
% time_red = [time_red,TI_NL.Solution.time];
% 
% 
% tmin = time_red(end);
% 
% ii = ii+1;
% 
% end
% 
% dyn.nlin_red_c.disp = decodeDofsNodes(u_full_red,nNodes,nDofPerNode); % (node, dof of node, tsamp)
% dyn.nlin_red_c.time = time_red;

%%
% Initial condition: equilibrium
u1_F = zeros(nDofsC,1);
v1_F = zeros(nDofsC,1);
a1_F = zeros(nDofsC,1);

tmin = 0;
tmax = tint;
u_full = [];
time_red = [];
u_full_red = [];
ud_full_red = [];

t0_ii = tmin;

pchange = [];
tchange = [];

while t0_ii < tmax

%choose which rom to use
p0_ii = k_t(t0_ii);
pchange = [pchange,p0_ii];
tchange = [tchange,t0_ii];

for jj = 1:size(k_validity_range,2)
    
    if (p0_ii <= k_validity_range(2,jj)) && (p0_ii >= k_validity_range(1,jj))
        
        kmin_ii = k_validity_range(1,jj);
        kmax_ii = k_validity_range(2,jj);
        break
        
    end
    
end
  
% q0_ii = ROMS{jj}.V'*u1_F;   
% qd0_ii = ROMS{jj}.V'*v1_F; 
% qdd0_ii = ROMS{jj}.V'*a1_F; 

pseudo_inv = pinv(ROMS{jj}.V);
q0_ii = pseudo_inv*u1_F;   
qd0_ii = pseudo_inv*v1_F; 
qdd0_ii = pseudo_inv*a1_F; 

%condition to terminate integration
term_cond = @(q,qd,t) exit_integration(q,qd,t,k_t,kmin_ii,kmax_ii);

% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark_mod('timestep',h,'alpha',0.005);

% nonlinear red Residual evaluation function handle;                  
residual = @(q,qd,qdd,t) residual_beam_spring_ROM(q,qd,qdd,t,BeamAssembly,F_ext,k_t,n_dof_spring,ROMS{jj}.V);

% integrate equations with Newmark scheme
TI_NL.Integrate(q0_ii,qd0_ii,qdd0_ii,t0_ii,tmax,residual,term_cond);

% store solution in more suitable array (only displacements)
u_red = TI_NL.Solution.q;
ud_red = TI_NL.Solution.qd;

u_full_ii = get_full_solution(u_red,ROMS{jj}.V);
ud_full_ii = get_full_solution(ud_red,ROMS{jj}.V);

%save final displacements and velocities
u1_F = u_full_ii(:,end);
v1_F = ROMS{jj}.V*(TI_NL.Solution.qd(:,end));
a1_F = ROMS{jj}.V*TI_NL.Solution.qdd_end;
%accelerations ?? 

%save full disp and time vector
u_full_ii = BeamAssembly.unconstrain_vector(u_full_ii); 
u_full_red = [u_full_red,u_full_ii];

ud_full_ii = BeamAssembly.unconstrain_vector(ud_full_ii); 
ud_full_red = [ud_full_red,ud_full_ii];

time_red = [time_red,TI_NL.Solution.time];


t0_ii = time_red(end);


end

dyn.nlin_red_c.disp = decodeDofsNodes(u_full_red,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn.nlin_red_c.vel = decodeDofsNodes(ud_full_red,nNodes,nDofPerNode); % (node, dof of node, tsamp)

dyn.nlin_red_c.time = time_red;

%% plot figure

%displacements
plot_dof_location = 0.9;%0.9%0.75;%2/3; %percentage of length of beam
node2plot = find_node(plot_dof_location*l,0,[],Nodes); % node where to put the force
plot_dof = get_index(node2plot, BeamAssembly.Mesh.nDOFPerNode );

figure;
%subplot 311; 
hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,1,:)),'-');
plot(dyn.nlin_red.time,squeeze(dyn.nlin_red.disp(node2plot,1,:)),'-');
plot(dyn.nlin_red_c.time,squeeze(dyn.nlin_red_c.disp(node2plot,1,:)),'-');
plot(tchange,zeros(length(tchange),1),'x','markersize',5,'color','r');
xlabel('t'); ylabel('u'); grid on; axis tight; 
title('dispx')

figure
%subplot 312; 
hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,2,:)),'-');
plot(dyn.nlin_red.time,squeeze(dyn.nlin_red.disp(node2plot,2,:)),'-');
plot(dyn.nlin_red_c.time,squeeze(dyn.nlin_red_c.disp(node2plot,2,:)),'-');
plot(tchange,zeros(length(tchange),1),'x','markersize',5,'color','r');
xlabel('t'); ylabel('v'); grid on; axis tight; 
title('dispy')

figure
%subplot 313; 
hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,3,:)),'-');
plot(dyn.nlin_red.time,squeeze(dyn.nlin_red.disp(node2plot,3,:)),'-');
plot(dyn.nlin_red_c.time,squeeze(dyn.nlin_red_c.disp(node2plot,3,:)),'-');
plot(tchange,zeros(length(tchange),1),'x','markersize',5,'color','r');
xlabel('t'); ylabel('w'); grid on; axis tight;
title('rot')

%velocities
figure;
subplot 311; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.vel(node2plot,1,:)),'-');
plot(dyn.nlin_red.time,squeeze(dyn.nlin_red.vel(node2plot,1,:)),'-');
plot(dyn.nlin_red_c.time,squeeze(dyn.nlin_red_c.vel(node2plot,1,:)),'-');
plot(tchange,zeros(length(tchange),1),'x','markersize',5,'color','r');
xlabel('t'); ylabel('u'); grid on; axis tight; 

subplot 312; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.vel(node2plot,2,:)),'-');
plot(dyn.nlin_red.time,squeeze(dyn.nlin_red.vel(node2plot,2,:)),'-');
plot(dyn.nlin_red_c.time,squeeze(dyn.nlin_red_c.vel(node2plot,2,:)),'-');
plot(tchange,zeros(length(tchange),1),'x','markersize',5,'color','r');
xlabel('t'); ylabel('v'); grid on; axis tight; 

subplot 313; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.vel(node2plot,3,:)),'-');
plot(dyn.nlin_red.time,squeeze(dyn.nlin_red.vel(node2plot,3,:)),'-');
plot(dyn.nlin_red_c.time,squeeze(dyn.nlin_red_c.vel(node2plot,3,:)),'-');
plot(tchange,zeros(length(tchange),1),'x','markersize',5,'color','r');
xlabel('t'); ylabel('w'); grid on; axis tight;

%% define error measure
% a = error_perc(dyn.nlin,dyn.nlin_red_c);
% 
% function err = error_perc(ref_sol,approx_sol,dofs,nodes)
% 
% 
% if nargin < 3
%     dofs = 1:1:size(approx_sol.disp,2);
%     nodes = 1:1:size(approx_sol.disp,3);
% end
% 
% n_time_samples = length(approx_sol.time);
% 
% ref_sol_disp = zeros(length(dofs)*length(nodes),n_time_samples);
% approx_sol_disp = zeros(length(dofs)*length(nodes),n_time_samples);
% 
% for ii = 1:length(nodes)
%     
%     ref_disp_nodeii = ref_sol.disp(nodes(ii),dofs,:);
%     ref_sol_disp()
%     
% end
% approx_disp = interp1(approx_sol.time,(approx_sol.disp).',ref_sol.time);
% 
% 
% end


%% functions


function [K,F] = Ktg_force(beamAssembly,k_spring,n_dof_spring, u)

[K, F] = beamAssembly.tangent_stiffness_and_force(u);
K(n_dof_spring,n_dof_spring) = +k_spring + K(n_dof_spring,n_dof_spring);
F(n_dof_spring) = +k_spring*u(n_dof_spring) + F(n_dof_spring);

end



function [ r, drdqdd,drdqd,drdq, c0] = residual_beam_spring( q, qd, qdd, t, beamAssembly, Fext, k_spring ,n_dof_spring)

k_spring_t = k_spring(t);

M = beamAssembly.DATA.M;
C = beamAssembly.DATA.C;

u = beamAssembly.unconstrain_vector(q);
[K, F] = Ktg_force(beamAssembly,k_spring_t,n_dof_spring, u); 
C_red = beamAssembly.constrain_matrix(C);
K_red = beamAssembly.constrain_matrix(K);
M_red = beamAssembly.constrain_matrix(M);
F_elastic = beamAssembly.constrain_vector(F);
F_external =  beamAssembly.constrain_vector(Fext(t)); 

% Residual is computed according to the formula above:
F_inertial = M_red * qdd;
F_damping = C_red * qd;
r = F_inertial + F_damping + F_elastic - F_external ;
drdqdd = M_red;
drdqd = C_red;
drdq = K_red;

c0 = norm(F_inertial) + norm(F_damping) + norm(F_elastic) + norm(F_external);

end



function [ROMS,k_ROMs,validity_range] = construct_ROM(k_samples,Assembly,nDofsF,n_dof_spring,nVMs,logic_MD)

nROMs = length(k_samples);
ROMS = cell(nROMs,1);

M = Assembly.mass_matrix();
M_C = Assembly.constrain_matrix(M);

u_eq = zeros(nDofsF,1);

for ii = 1:nROMs
   
    [K,~] = Ktg_force(Assembly,k_samples(ii),n_dof_spring, u_eq);
    
    K_C = Assembly.constrain_matrix(K); %tangent stiffness matrix for u = 0
    
    [VM,omega2] = eigs(K_C,M_C,nVMs,'SM'); 
    %omega = sqrt(diag(omega2));
    %normalize eigenvectors
    for jj = 1:nVMs
        VM(:,jj) = VM(:,jj)/vecnorm(VM(:,jj));
    end
    
    
    %compute here the modal derivatives
    if logic_MD == 1
        
        %compute Modal Derivatives ii
        [MD, names] = modal_derivatives(Assembly, VM, K_C);
        
        for jj = 1:size(MD,2)
            MD(:,jj) = MD(:,jj)./vecnorm(MD(:,jj)); %normalize eigenvectors (they are unitary vectors)
        end
        
    end
    
    if logic_MD == 1
        ROMS{ii}.V = orth([VM,MD]);
    else
        ROMS{ii}.V = VM;
    end
    
    
end 
    

k_ROMs = k_samples; 


aa = [k_ROMs(1),k_ROMs,k_ROMs(end)];
bb = (aa(1:end-1) + aa(2:end))/2;

validity_range = zeros(2,length(bb)-1);

for ii = 1:size(validity_range,2)
    validity_range(:,ii) = [bb(ii),bb(ii+1)].';
end
    
end


function [ r, drdqdd,drdqd,drdq, c0] = residual_beam_spring_ROM( q, qd, qdd, t, beamAssembly, Fext, k_spring ,n_dof_spring,V)

k_spring_t = k_spring(t);

M = beamAssembly.DATA.M;
C = beamAssembly.DATA.C;


u = beamAssembly.unconstrain_vector(V*q);
[K, F] = Ktg_force(beamAssembly,k_spring_t,n_dof_spring, u); 

C_red = V.'*beamAssembly.constrain_matrix(C)*V;
K_red = V.'*beamAssembly.constrain_matrix(K)*V;
M_red = V.'*beamAssembly.constrain_matrix(M)*V;

F_elastic = V.'*beamAssembly.constrain_vector(F);
F_external =  V.'*beamAssembly.constrain_vector(Fext(t)); 

% Residual is computed according to the formula above:
F_inertial = M_red * qdd;
F_damping = C_red * qd;
r = F_inertial + F_damping + F_elastic - F_external ;
drdqdd = M_red;
drdqd = C_red;
drdq = K_red;

c0 = norm(F_inertial) + norm(F_damping) + norm(F_elastic) + norm(F_external);

end

function x_full = get_full_solution(disp_red,V)

n_samples_t = size(disp_red,2);
dim_full = size(V,1);
x_full = zeros(dim_full,n_samples_t);

for ii = 1:n_samples_t
   
    x_full(:,ii) = V*disp_red(:,ii);
    
end

end


function logic = exit_integration(q,qd,t,k_t,kmin,kmax)

if k_t(t) < kmin || k_t(t) > kmax
   logic = 1;
else
   logic = 0;
end

end