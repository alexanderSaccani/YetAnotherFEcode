% author: ALEXANDER SACCANI - PHD CANDIDATE, ETH ZURICH, 
% created: 1.12.2022
% last modified: 30.01.2023

% Straight beam clamped at both ends

clearvars ; close all; clc
%% Define colors for plots

%% Model construction

%GEOMETRY__________________________________________________________________
l = 0.1; 
h = 1e-3;
b = 1e-2;

height_midspan = 2e-3;
w = height_midspan/l; % this parameter varies from 0 to 0.5 (straight beam to half circumference)

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


%create mesh from nodes
Nodes = [nodes_x, nodes_y];
arcMesh = Mesh(Nodes); 

%define elements 
Elements = [1:nNodes-1;2:nNodes].'; % ROW -> # element, 
                                    % COLUMN -> # of node connected
myElementConstructor = @()ThermalBeamElement(b, h, myThermalBeamMaterial); %type of element

%create element table from elements
arcMesh.create_elements_table(Elements,myElementConstructor);


%BOUNDARY CONDITIONS_______________________________________________________
arcMesh.set_essential_boundary_condition([1 arcMesh.nNodes],[1 2 3],0) % Doubly clamped beam
%BeamMesh.set_essential_boundary_condition(1,[1 2 3],0) % clamped beam  (node,[constrained dofs],imposed disp.)

BeamAssembly = Assembly(arcMesh);

nDofsF = BeamAssembly.Mesh.EBC.nDOFs; % #dofs free (without bc)
nDofsC = nDofsF - size(BeamAssembly.Mesh.EBC.constrainedDOFs,1); % #dofs constrained (with bc)


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
    num2str(omega(mode2plot(ii))/(2*pi)) ' Hz'])
    
    subplot 412
    plot(nodes_x,v_ii,'-'); title('tranverse disp.'); xlabel('x')
    subplot 413
    plot(nodes_x,u_ii,'-'); title('longitudinal disp.'); xlabel('x')
    subplot 414
    plot(nodes_x,r_ii,'-'); title('rotation'); xlabel('x')
    
end


%% Define loads 
% Nodal force
F_c = zeros(nDofsF,1);  %concentrated load

% transverse load 
fext_node = find_node(l/2,0,[],Nodes); % node where to put the force
fext_dofs = get_index(fext_node, nDofPerNode);
F_c(fext_dofs(2)) = 3;%30; 

%uniform body force
Pressure_v = 0; 
F_v_uniform = Pressure_v*BeamAssembly.uniform_body_force();
F_v_uniform(1:3:end) = F_v_uniform(1:3:end)*0; %set to 0 the forces corresponding to the axial dofs. (no axial pressure)

%define total external load
F = F_c + F_v_uniform;

%% Static Analysis 
% Define static thermal sin^2 distribution
T_ampl = 30;
p = l/4; %width of thermal pulse

xc_st = l/2; %center of thermal pulse
x0_st = xc_st - p/2; %left extreme of pulse
T_st = T_ampl*(sin(pi*(nodes_x-x0_st)/p)).^2.*(heaviside(nodes_x-x0_st) - ...
    heaviside((nodes_x-x0_st)-p)); %static temp profile

% linear static solution___________________________________________________
[K_hot,gth] = BeamAssembly.tangent_stiffness_and_force(zeros(nDofsF,1),T_st);
u_lin = BeamAssembly.solve_system(K_hot, F-gth);   %subtract the internal force generated  by the temperature
static.lin = reshape(u_lin,3,[]).';	% Linear response

% nonlinear static solution________________________________________________
BeamAssembly.DATA.K = K_hot;    % this value of K is used only for first guess solution 
[~,u] = static_equilibrium_thermal(BeamAssembly, F, T_st,'display', 'iter-detailed');
static.nlin = reshape(u,nDofPerNode,[]).';     % Nonlinear response

%plot static displacements_________________________________________________
scl = 50;   %scaling factor

figure('units','normalized','position',[.1 .1 .8 .8])
subplot 411; 
plot(nodes_x,T_st,'linewidth',2); xlabel('x [m]'); ylabel('T [K]'); 
title('Static analysis. T distribution')

subplot 412; 
plot(nodes_x,static.nlin(:,1),'color','r','linewidth',2); hold on;
plot(nodes_x,static.lin(:,1),'--', 'color','cyan','linewidth',2)
xlabel('x [m]'); ylabel('u [m]'); title('long. displ');legend('nonlinear','linear','linewidth',2)

subplot 413; 
plot(nodes_x,static.nlin(:,2),'color','r','linewidth',2); hold on;
plot(nodes_x,static.lin(:,2),'--', 'color','cyan','linewidth',2);
xlabel('x [m]'); ylabel('v [m]'); title('tranv. displ'); legend('nonlinear','linear','linewidth',2)

subplot 414; 
plot(nodes_x,nodes_y,'.-','markersize',10,'color','b','linewidth',2); hold on;
plot(nodes_x + scl*static.nlin(:,1), nodes_y + scl*static.nlin(:,2), '.-', 'markersize',10,'color','r','linewidth',2)
plot(nodes_x + scl*static.lin(:,1), nodes_y + scl*static.lin(:,2), '--', 'color','cyan','linewidth',2)
title('static deformation shape'); xlabel('x [m]'); legend('undeformed','nonlinear','linear','linewidth',2)

%% static Buckling analysis
F0 = F; %load configuration corresponding to lambda = 1 (load multiplier)

lam_high = 12; %load multiplier max abs value (symmetric analysis for positive and negative values of lambda from the origin

%T samples 
T_samples = [0,50,100,160];

%parameters for continuation
Sopt.parametrization = 'arc_length';
Sopt.predictor = 'secant';
Sopt.reversaltolerance = 1;
ds = {0.03,0.04,0.03,0.01};

%initialize output
X_buckl_an = cell(length(T_samples),1);

for ii = 1:length(T_samples)
    
    T_ampl_ii = T_samples(ii);
    
    p = l/4; %width of thermal pulse
    xc_st = l/2; %center of thermal pulse
    x0_st = xc_st - p/2; %left extreme of pulse
    
    T = T_ampl_ii*(sin(pi*(nodes_x-x0_st)/p)).^2.*(heaviside(nodes_x-x0_st) - ...
        heaviside((nodes_x-x0_st)-p)); %static temp profile


    [K_hot,gth] = BeamAssembly.tangent_stiffness_and_force(zeros(nDofsF,1),T);
    u_lin_0 = BeamAssembly.solve_system(K_hot, 0-gth);   %subtract the internal force generated  by the temperature

    u_guess = BeamAssembly.constrain_vector(u_lin_0);

    % ds = norm(u_lin_0);
    
    ds_ii = ds{ii};
    
    fun_residual = @(X) residual_buckling(X,BeamAssembly,F0,T);
    [X1,Solinfo,Sol] = solve_and_continue(u_guess,fun_residual,0,lam_high,ds_ii,Sopt);
    [X2,Solinfo,Sol] = solve_and_continue(u_guess,fun_residual,0,-lam_high,ds_ii,Sopt);
    
    X_buckl_an{ii} = [flip(X2,2),X1];

       
end

% plot results of buckling analysis
plot_dof_location = 0.5; %percentage of length of beam
node2plot = find_node(plot_dof_location*l,0,[],Nodes); % node to plot
dof2plot = get_index(node2plot, nDofPerNode);
dof2plot = dof2plot(2); %plot vertical displacement

color_list = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE',...
    '#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};

figure('units','normalized','position',[.3 .3 .4 .4]);
title(['buckling analysis: tranversal displacement at ', num2str(plot_dof_location*100),'% of beam span']);
leg = cell(length(T_samples),1); %legend initialize
for ii = 1:length(T_samples)
    
    hold on;
    X_ii = X_buckl_an{ii};
    delta = X_ii(dof2plot,:);
    lam = X_ii(end,:);
    plot(delta,lam,'color',color_list{ii},'linewidth',2);
    
    xlabel('displacement [m]')
    ylabel('\lambda')
    
    leg{ii} = ['T = ',num2str(T_samples(ii))];
    
end
legend(leg);
grid on;
%% Dynamic response using Implicit Newmark
%force and thermal load____________________________________________________
%define external forcing
om_fext = (omega(1)+omega(2))/2;%4e3; %omega(1)/10;% + (omega(2) - omega(1))/12;
T_fext = 2*pi/om_fext; % period of ex. forcing

%construct forcing vector
F_ext = @(t) F*sin(om_fext*t); 

% % 1)thermal load tranversing the beam
% %define thermal load 
% T_ampl = 30;
% p = l/4; %width of thermal pulse
% A = l+p;  %amplitude of oscillation of pulse
% xci = -p/2; %center location at initial time
% eps = 0.01; %om_fex/om_th (th = thermal forcing)
% om_th = eps*om_fext;  % angular frequency of thermal mode oscillation in time
% 
% T_th = 2*pi/om_th; % period of temp. oscillations
% 
% x0 = @(xc) xc - p/2; %left extreme of pulse
% T_dyn_p = @(xc) T_ampl*(sin(pi*(nodes_x-x0(xc))/p)).^2.*(heaviside(nodes_x-x0(xc)) - ...
%     heaviside((nodes_x-x0(xc))-p)); %define the temperature as a function of xc only (T profile parametrized with xc
% 
% xc_dyn = @(t) xci + A*sin(om_th*t); %define how center of pulse xc variates in time
% 
% T_dyn_t = @(t) T_dyn_p(xc_dyn(t)); %evaluate the temperature profile in time

% 2)thermal pulse at the center ramping up
T_ampl = 160;
p = l/2; %width of thermal pulse
xc = l/2; %center location at initial time
eps = 0.05; %om_fex/om_th (th = thermal forcing)
om_th = eps*om_fext;  % angular frequency of thermal mode oscillation in time

T_th = 2*pi/om_th; % period of temp. oscillations

x0 = xc - p/2;
T_dyn_p = @(T_ampl) T_ampl*(sin(pi*(nodes_x-x0)/p)).^2.*(heaviside(nodes_x-x0) - ...
    heaviside((nodes_x-x0)-p));

T_ampl_t = @(t) T_ampl/T_th*4*t - (heaviside(t-T_th/4))*T_ampl/T_th*4*(t-T_th/4);

T_dyn_t = @(t) T_dyn_p(T_ampl_t(t));



% %plot temperature profile (animation)
% time = linspace(0,T_th/4+2.5*T_th/4,100);
% figure
% for ii = 1:length(time)
%     plot(nodes_x,T_dyn_t(time(ii)));
%     axis([0,l,0,T_ampl]);
%     title(['T distribution, t = ',num2str(time(ii)),'s']);
%     pause(0.03);
% end


%add step term to forcing
t_snap = 1.5*T_th/4;
delta_t = 7*T_th/4;
mult = 25; 

F_ext = @(t)  (- mult*F*heaviside(t-t_snap) + mult*F*heaviside(t-t_snap-delta_t));
time = linspace(0,tint,100);
figure
plot(time,F_ext(time));

%% Thermal VMs analysis (just to visualize them)
% p_sampl = [-p/2;linspace(p/2,l-p/2,3)']; %sample locations of center thermal pulse

p_sampl = linspace(0,T_ampl,5)'; %sample locations of center thermal pulse

T_sampl = zeros(nNodes,length(p_sampl)); % T nodal distribution samples

%generate corresponding temperature profiles
for ii = 1:length(p_sampl)
    T_sampl(:,ii) = T_dyn_p(p_sampl(ii));
end

VMs_at_eq = 1; %do you want to compute VMs at thermal equilibrium? (set to 1 if yes, otherwise set to 0)
[VMs,static_eq,omega] = VM_thermal(BeamAssembly,T_sampl,VMs_at_eq,2);

%plot VMs__________________________________________________________________
T_sampl_2_plot = 1:length(p_sampl);
mode2plot = [1,2];

%T distribution 
figure('units','normalized','position',[.3 .3 .4 .4]); hold on;
title('T distribution samples'); xlabel('x [m]'); ylabel('T [K]');
for jj = 1:length(T_sampl_2_plot)
    plot(nodes_x,T_sampl(:,T_sampl_2_plot(jj)))
end

%VMs
for ii = 1:length(mode2plot)
    
    VM2pl = mode2plot(ii);
    figure('units','normalized','position',[.3 .3 .4 .4]); hold on;
    title(['VM ',num2str(mode2plot(ii))]); xlabel('x [m]');
    
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
   title('Static thermal equilibrium '); xlabel('x [m]');
   ueq_ii = reshape(static_eq(:,ii),nDofPerNode,[]).';
   subplot 311; hold on;
   plot(nodes_x,ueq_ii(:,1)); 
   xlabel('x [m]'); ylabel('u [m]');
end
for ii = 1:length(T_sampl_2_plot)
   ueq_ii = reshape(static_eq(:,ii),nDofPerNode,[]).';
   subplot 312; hold on;
   plot(nodes_x,ueq_ii(:,2)); 
   xlabel('x [m]'); ylabel('v [m]');
end
scl = 50; %scaling
for ii = 1:length(T_sampl_2_plot)
   ueq_ii = reshape(static_eq(:,ii),nDofPerNode,[]).';
   subplot 313; hold on;
   plot(nodes_x + scl*ueq_ii(:,1), nodes_y + scl*ueq_ii(:,2),'*-'); 
   xlabel('x [m]'); ylabel('deformation')
end

%% integration of EOMs

% Integrate linearized model_______________________________________________
% settings for integration
% tint = T_th/4; % integration interval

tint = 12*T_th/4;%+5*T_th/4; % integration interval

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

% linear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_thermal_linear(q,qd,qdd,t,...
    BeamAssembly, F_ext, T_dyn_t);

% integrate equations with Newmark scheme
TI_LIN.Integrate(q0,qd0,qdd0,tint,residual);

% store solution in more suitable array (only displacements)
u_dyn = BeamAssembly.unconstrain_vector(TI_LIN.Solution.q); 
dyn.lin.disp = decodeDofsNodes(u_dyn,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn.lin.time = TI_LIN.Solution.time;

% Integrate nonlinear EOM__________________________________________________
% settings for integration

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
TI_NL.Integrate(q0,qd0,qdd0,tint,residual);

% store solution in more suitable array (only displacements)
u_dyn = BeamAssembly.unconstrain_vector(TI_NL.Solution.q); 
dyn.nlin.disp = decodeDofsNodes(u_dyn,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn.nlin.time = TI_NL.Solution.time;

% Quasistatic Solution_____________________________________________________
dyn.quasistatic.time = linspace(0,tint,1e2); %define time vector
u_quasist = quasistatic_solution( BeamAssembly, T_dyn_t, dyn.quasistatic.time ); %compute quasistatic solution (F=0)
dyn.quasistatic.disp = decodeDofsNodes(u_quasist,nNodes,nDofPerNode); % (node, dof of node, tsamp)

% plot results of nonlinear dynamic analysis_______________________________
plot_dof_location = 0.3; %percentage of length of beam
node2plot = find_node(plot_dof_location*l,0,[],Nodes); % node to plot

figure('units','normalized','position',[.1 .1 .8 .8]);
subplot 311; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,1,:)),'-','color','r','linewidth',2);
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,1,:)),'--','color','cyan','linewidth',2);
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,1,:)),'--','color','green','linewidth',2);
xlabel('t [s]'); ylabel('u [m]'); grid on; axis tight;  legend('nonlinear','linear','quasistatic','linewidth',2);

subplot 312; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,2,:)),'-','color','r','linewidth',2);
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,2,:)),'--','color','cyan','linewidth',2);
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,2,:)),'--','color','green','linewidth',2);
xlabel('t [s]'); ylabel('v [m]'); grid on; axis tight; legend('nonlinear','linear','quasistatic','linewidth',2);

subplot 313; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,3,:)),'-','color','r','linewidth',2);
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,3,:)),'--','color','cyan','linewidth',2);
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,3,:)),'--','color','green','linewidth',2);
xlabel('t [s]'); ylabel('w [-]'); grid on; axis tight;  legend('nonlinear','linear','quasistatic');
%%
 %plot results in time (animation)______________________________________
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
%    plot(nodes_x,T_dyn_t(time_tmp(ii)));
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
%

% scl = 20; %scale factor
% def_ax_bounds = [0,l,1.1*(min(min((scl*v_tmp)))-min(nodes_y)),...
%     1.1*(max(max(scl*v_tmp))+max(nodes_y))];
% 
% figure('units','normalized','position',[.1 .1 .8 .8]);
% for ii = 1:length(time_tmp)
%     
%     time_ii = time_tmp(ii);
%     
%     subplot 211
%     plot(nodes_x,T_dyn_t(time_tmp(ii)),'linewidth',1.5,'color','r');
%     axis([0,l,0,T_ampl]);
%     xlabel('x [m]'); ylabel('T [K]');
%     title('temperature profile')
% 
%     subplot 212
%     plot(nodes_x + scl*u_tmp(:,ii), nodes_y + scl*v_tmp(:,ii), 'marker', '.','markersize',10,'color','b','linewidth',1.5 );
%     hold on;
%     
%     plot(nodes_x, nodes_y, 'marker', '.','color','k' ,'linewidth',1.5,'markersize',10);
%     axis(def_ax_bounds)
%     xlabel('x [m]'); ylabel('y [m]');
%     title(['arc vibration, time ', num2str(round(time_ii,4)),' s, scale factor: ',num2str(scl) ])
%     legend('deformed','undeformed')
%     
%     pause(0.001);
%     
%     hold off;
% end



%% Reduced Order Model Construction (constant basis)

T_ampl_sampl = 0;

T_samples = T_dyn_p(T_ampl_sampl);


%construct model___________________________________________________________
% define how many VMs to use in modal reduction (also the corresponding
% modal derivatives are added)
number_VMs = 3;

logic_MD = 1;

out = VMMD_thermal(BeamAssembly,T_samples,1,number_VMs,logic_MD);

VMs = BeamAssembly.constrain_vector(out.VMs{1});
MDs = BeamAssembly.constrain_vector(out.MDs{1});

V = [VMs,MDs];

V = orth(V);


%% Integrate nonlinear ROM (constant basis)
% settings for integration
h = T_fext/50; % time step for integration

% Initial condition: equilibrium
n_dof_r = size(V,2);

q0 = zeros(n_dof_r,1);
qd0 = zeros(n_dof_r,1);
qdd0 = zeros(n_dof_r,1);

% Instantiate object for nonlinear time integration
TI_NL_r = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear Residual evaluation function handle
residual_r = @(q,qd,qdd,t)residual_ROM(q,qd,qdd,t,BeamAssembly,V,F_ext,T_dyn_t);

% integrate equations with Newmark scheme
TI_NL_r.Integrate(q0,qd0,qdd0,tint,residual_r);

% postprocess output 
u_dyn_r = V*TI_NL_r.Solution.q;

u_dyn_r = BeamAssembly.unconstrain_vector(u_dyn_r); 
dyn_r.nlin.disp = decodeDofsNodes(u_dyn_r,nNodes,nDofPerNode); % (node, dof of node, tsamp)
dyn_r.nlin.time = TI_NL_r.Solution.time;


%% Plot comparison
% plot results of nonlinear dynamic analysis_______________________________
plot_dof_location = 0.4; %percentage of beam length
node2plot = find_node(plot_dof_location*l,0,[],Nodes); % node where to put the force
plot_dof = get_index(node2plot, BeamAssembly.Mesh.nDOFPerNode );

figure('units','normalized','position',[.1 .1 .8 .8]); 
subplot 311; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,1,:)),'-','color','r','linewidth',2);
plot(dyn_r.nlin.time,squeeze(dyn_r.nlin.disp(node2plot,1,:)),'--','color','b','linewidth',2);
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,1,:)),'--','color','cyan','linewidth',2);
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,1,:)),'--','color','g','linewidth',2);
xlabel('t [s]'); ylabel('u [m]'); grid on; axis tight; legend('full','ROM','lin','quasistatic');
title(['comparison of ROM and HFM solutions, dof at ',num2str(plot_dof_location*100), '% of beam length'])

subplot 312; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,2,:)),'-','color','r','linewidth',2);
plot(dyn_r.nlin.time,squeeze(dyn_r.nlin.disp(node2plot,2,:)),'--','color','b','linewidth',2);
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,2,:)),'--','color','cyan','linewidth',2);
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,2,:)),'--','color','g','linewidth',2);
xlabel('t [s]'); ylabel('v [m]'); grid on; axis tight; legend('full','ROM','lin','quasistatic');

subplot 313; hold on;
plot(dyn.nlin.time,squeeze(dyn.nlin.disp(node2plot,3,:)),'-','color','r','linewidth',2);
plot(dyn_r.nlin.time,squeeze(dyn_r.nlin.disp(node2plot,3,:)),'--','color','b','linewidth',2);
plot(dyn.lin.time,squeeze(dyn.lin.disp(node2plot,3,:)),'--','color','cyan','linewidth',2);
plot(dyn.quasistatic.time,squeeze(dyn.quasistatic.disp(node2plot,3,:)),'--','color','g','linewidth',2);
xlabel('t [s]'); ylabel('w [-]'); grid on; axis tight; legend('full','ROM','lin','quasistatic');

%% residual function for ROM

function [ r, drdqdd,drdqd,drdq, c0] = residual_ROM( q, qd, qdd, t, Assembly, V, Fext, T_fn_t)


%compute T((t))
T = T_fn_t(t);


M = Assembly.constrain_matrix(Assembly.DATA.M); %can be more efficient by providing the constrain matrix (then pay attention to ROMs)
C = Assembly.constrain_matrix(Assembly.DATA.C);

Mr = V.'*M*V; %constrain matrix
Cr = V.'*C*V; %constrain matrix

u = Assembly.unconstrain_vector(V*q); %get full coordinates

[Ktg, Fi_full] = Assembly.tangent_stiffness_and_force(u,T); %get full internal force and tangent stiffness
F_int_r = V.'*Assembly.constrain_vector(Fi_full); %back to reduced space
Ktgr = V.'*Assembly.constrain_matrix(Ktg)*V;

%compute residual
F_inertial = Mr*qdd;
F_damping = Cr*qd;
F_ext_r = V.'*Assembly.constrain_vector(Fext(t));

r = F_inertial + F_damping + F_int_r - F_ext_r;

%Jacobians
drdqdd = Mr;
drdqd = Cr;
drdq = Ktgr;
%% 
% We use the following measure to comapre the norm of the residual $\mathbf{r}$
% 
% $$\texttt{c0} = \|\mathbf{M_V}\ddot{\mathbf{q}}\| + \|\mathbf{C_V}\dot{\mathbf{q}}\| 
% + \|\mathbf{F_V}(\mathbf{q})\| + \|\mathbf{V}^{\top}\mathbf{F}_{ext}(t)\|$$
c0 = norm(F_inertial) + norm(F_damping) + norm(F_int_r) + norm(F_ext_r);
end














