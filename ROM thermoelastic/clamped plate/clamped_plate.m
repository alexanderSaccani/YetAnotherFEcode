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
height_midspan = 0.0000001; %height of midspan 

%Mesh______________________________________________________________________
nx = 15; % # elements along x
ny = 8; % # elements along y

eltype = 'QUAD8'; %'QUAD4'
myElementConstructor = @()Quad8Shell(thickness, myMaterial,2);
[nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,eltype);

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

nNodes = myMesh.nNodes;

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
PlotMesh([nodes,zNodes],elements);

%Boundary conditions: all sides are clamped________________________________
myMesh.set_essential_boundary_condition([nset{1}],1:5,0)
myMesh.set_essential_boundary_condition([nset{2}],1:5,0)
myMesh.set_essential_boundary_condition([nset{3}],1:5,0)
myMesh.set_essential_boundary_condition([nset{4}],1:5,0)

%Assembly__________________________________________________________________
Assembly = Assembly(myMesh);
M = Assembly.mass_matrix();
u0 = zeros( myMesh.nDOFs, 1);

%Compute and store Mass matrix, K cold_____________________________________
%T field (zeros)
T = zeros(nNodes,1);
gradT = zeros(nNodes,1);
[Kc,~] = Assembly.tangent_stiffness_and_force(u0,T,gradT);

% store matrices
Assembly.DATA.Kc = Kc; %Kcold
Assembly.DATA.M = M; %mass matrix


%% Vibration modes (cold)
nVMs = 5; % first n_VMs modes with lowest frequency calculated 

Kconstr = Assembly.constrain_matrix(Assembly.DATA.Kc);
Mconstr = Assembly.constrain_matrix(Assembly.DATA.M);

[Vcold,om] = eigs(Kconstr, Mconstr, nVMs, 'SM');
[fcold,ind] = sort(sqrt(diag(om))/2/pi);

Vcold = Vcold(:,ind);
for ii = 1:nVMs
    Vcold(:,ii) = Vcold(:,ii)/max(sqrt(sum(Vcold(:,ii).^2,2)));
end
Vcold = Assembly.unconstrain_vector(Vcold);

%plot VMs__________________
VMs = cell(nVMs,1);
for ii = 1:nVMs  
    VMs{ii,1} = reshape(Vcold(:,ii),5,[]).';    
end

for ii = 1:nVMs
    figure; hold on;
    VM2plot = VMs{ii};
    PlotFieldonDeformedMesh([nodes,zNodes],elements, [VM2plot(:,1),...
        VM2plot(:,2), VM2plot(:,3)], 'factor',max(nodes(:,2)));
    title(['VM ',num2str(ii),', f: ', num2str(fcold(ii)),' Hz']);
end

%% Buckling Analysis for different temperatures

BUCKLING_ANALYSIS = 0; %set to 0 to use data stored in variable "filename"
filename = "bucklAnal.mat";

if BUCKLING_ANALYSIS == 1
    
    %define forcingF = zeros(myMesh.nDOFs,1);
    nf = find_node(Lx/2,Ly/2,[],nodes); % node where to put the force
    node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
    F(node_force_dofs(3)) = +1;

    F0 = F; %load configuration corresponding to lambda = 1 (load multiplier)

    lam_high = 20; %load multiplier max abs value (symmetric analysis for positive and negative values of lambda from the origin

    %parameters for continuation
    Sopt.parametrization = 'orthogonal';
    Sopt.predictor = 'tangent';
    Sopt.noconv_stepcor = 'red';
    Sopt.errmax = 6;
    Sopt.reversaltolerance = 1;
    ds = 0.001;
    Sopt.dsmin = ds/100;
    Sopt.dsmax = ds*1000;

    %T samples 
    T_samples = [0,5,10];
    gradT = zeros(myMesh.nNodes,1);

    %initialize output
    X_buckl_an = cell(length(T_samples),1);

    Solopt = optimset(optimset(@fsolve),'Display','off',...%'iter',...
            'Jacobian','on','MaxIter',50,'DerivativeCheck','off');

    for ii = 1:length(T_samples)

        T_ampl_ii = T_samples(ii);

        T = T_ampl_ii*ones(myMesh.nNodes,1); %static temp profile

        [K_hot,gth] = Assembly.tangent_stiffness_and_force(zeros(myMesh.nDOFs,1),T,gradT);
    %   u_lin_0 = Assembly.solve_system(K_hot, zeros(myMesh.nDOFs,1), gth);   %subtract the internal force generated  by the temperature

        u_lin_FAKE = Assembly.solve_system(Assembly.DATA.Kc, zeros(myMesh.nDOFs,1), gth); %compute the guess with cold tangent stiffness (otherwise it explodes)
    %   u_guess = Assembly.constrain_vector(u_lin_0);

        u_guess = Assembly.constrain_vector(u_lin_FAKE);
        % ds = norm(u_lin_0);

        fun_residual = @(X) residual_buckling(X,Assembly,F0,T,gradT);
        [X1,Solinfo,Sol] = solve_and_continue(u_guess,fun_residual,0,lam_high,ds,Sopt,[],Solopt);
        [X2,Solinfo,Sol] = solve_and_continue(u_guess,fun_residual,0,-lam_high,ds,Sopt,[],Solopt);

        X_buckl_an_ii = [flip(X2,2),X1];
        lam_ii = X_buckl_an_ii(end,:);
        X_ii = Assembly.unconstrain_vector(X_buckl_an_ii(1:end-1,:));

        X_buckl_an{ii} = [X_ii;lam_ii];


    end

%      save("bucklAnal.mat","T_samples","X_buckl_an","F0");

else
    
    load(filename)

end

% plot results of buckling analysis
plot_dof_location = [0.5,0.5]; %percentage of length of beam
node2plot = find_node(plot_dof_location(1)*Lx,plot_dof_location(2)*Ly,[],nodes); % node to plot
nDofPerNode = 5;
dof2plot = get_index(node2plot, nDofPerNode);
dof2plot = dof2plot(3); %plot vertical displacement

color_list = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE',...
    '#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};

figure('units','normalized','position',[.3 .3 .4 .4]);
title(['buckling analysis: tranversal displacement at [', ...
    num2str(plot_dof_location(1)*100),',',num2str(plot_dof_location(2)*100),']','% of panel span']);
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


%% Hot modes
T_ampls = [5,10];

n_sampl = length(T_ampls);
T_samples = zeros(myMesh.nNodes,n_sampl);
gradT_samples = zeros(myMesh.nNodes,n_sampl);

for ii = 1:n_sampl
    T_samples(:,ii) = ones(myMesh.nNodes,1)*T_ampls(ii);
end

nVMs = 5;
[VMsHot,static_eq,omega] = VM_thermal(Assembly,nVMs,T_samples,gradT_samples);
fhot = omega/2/pi;

%Plots_____________________________________________________________________
for jj = 1:length(T_ampls)
    
    VMs = cell(nVMs,1);
    for ii = 1:nVMs  
        VMs{ii,1} = reshape(VMsHot{jj}(:,ii),5,[]).';    
    end

    for ii = 1:nVMs
        figure; hold on;
        VM2plot = VMs{ii};
        PlotFieldonDeformedMesh([nodes,zNodes],elements, [VM2plot(:,1),...
            VM2plot(:,2), VM2plot(:,3)], 'factor',max(nodes(:,2)));
        title(['VM ',num2str(ii),' dT = ',num2str(T_ampls(jj)),'K, f: ', num2str(fhot(ii,jj)),' Hz']);
    end

end


return




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


return

%% Static analysis 
% T distribution
T = zeros(nNodes,1);
gradT = zeros(nNodes,1);

% Nodal force
F = zeros(myMesh.nDOFs,1);
nf = find_node(Lx/2,Ly/2,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
%F(node_force_dofs(2)) = +10e5;
F(node_force_dofs(3)) = +1;


%Temperature field
T = 0*ones(nNodes,1);
gradT = 0*randn(nNodes,1);

u_lin = Assembly.solve_system(K_c, F);
ULIN = reshape(u_lin,5,[]).';	% Linear response
[~,u] = static_eq_shell_thermal(Assembly, F, T, gradT, 'display', 'iter-detailed');
UNL = reshape(u,5,[]).';        % Nonlinear response


% 
% node_plot_dofs = get_index(nset{2}, myMesh.nDOFPerNode );
% dofs_u = node_plot_dofs(1:5:end);
% dofs_v = node_plot_dofs(2:5:end);
% dofs_w = node_plot_dofs(3:5:end);
% dofs_thx = node_plot_dofs(4:5:end);
% dofs_thy = node_plot_dofs(5:5:end);
% set2plot = 4;
% 
% x_plot = (nodes(nset{set2plot}))';
% 
% u = UNL(nset{set2plot},1);
% v = UNL(nset{set2plot},2);
% w = UNL(nset{set2plot},3);
% thx = UNL(nset{set2plot},4);
% thy = UNL(nset{set2plot},5);
% 
% ul = ULIN(nset{set2plot},1);
% vl = ULIN(nset{set2plot},2);
% wl = ULIN(nset{set2plot},3);
% thxl = ULIN(nset{set2plot},4);
% thyl = ULIN(nset{set2plot},5);


Ulinx = ULIN(:,1);
Uliny = ULIN(:,2);
Ulinz = ULIN(:,3);
% 
UNlinx = UNL(:,1);
UNliny = UNL(:,2);
UNlinz = UNL(:,3);


figure
PlotFieldonDeformedMesh([nodes,zNodes],elements, [Ulinx, Uliny, Ulinz], 'factor',12)
hold on;
PlotMesh([nodes,zNodes],elements);
grid on
title('linear static solution');
axis

figure
PlotFieldonDeformedMesh([nodes,zNodes],elements, [UNlinx, UNliny, UNlinz], 'factor',20)
hold on;
PlotMesh([nodes,zNodes],elements);
title('non linear static solution');
grid on
axis

%% Static analysis


%% load displacement curve

F = zeros(myMesh.nDOFs,1);
nf = find_node(Lx/2,Ly/2,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F(node_force_dofs(3)) = +10e4;


%Temperature field
T = 0*ones(nNodes,1);
gradT = 0*randn(nNodes,1);

F0 = F; %load configuration corresponding to lambda = 1 (load multiplier)

lam_high = 9; %load multiplier max abs value (symmetric analysis for positive and negative values of lambda from the origin

%parameters for continuation
Sopt.parametrization = 'arc_length';
Sopt.predictor = 'secant';
Sopt.reversaltolerance = 1;
ds = 0.1;


[K_hot,gth] = Assembly.tangent_stiffness_and_force(zeros(myMesh.nDOFs,1),T,gradT);
u_lin_0 = Assembly.solve_system(K_hot, 0-gth);   %subtract the internal force generated  by the temperature
u_guess = Assembly.constrain_vector(u_lin_0);

    % ds = norm(u_lin_0);
fun_residual = @(X) residual_buckling(X,Assembly,F0,T,gradT);
[X1,Solinfo,Sol] = solve_and_continue(u_guess,fun_residual,0,lam_high,ds,Sopt);
[X2,Solinfo,Sol] = solve_and_continue(u_guess,fun_residual,0,-lam_high,ds,Sopt);
    
X_buckl_an = [flip(X2,2),X1];


% plot results of buckling analysis
plot_dof_location = 0.5; %percentage of length of beam
node2plot = find_node(plot_dof_location*Lx/2,Ly/2,[],nodes); % node to plot
dof2plot = get_index(node2plot, 5);
dof2plot = dof2plot(3); %plot vertical displacement


delta = X_buckl_an(dof2plot,:);
lam = X_buckl_an(end,:);

figure
hold on;
plot(delta,lam,'color','b','linewidth',2);
xlabel('displacement [m]')
ylabel('\lambda')














return

figure 
plot(x_plot,u)
hold on
plot(x_plot,ul)
legend('nlin','lin')

figure 
plot(x_plot,w)
hold on
plot(x_plot,wl)
legend('nlin','lin')


figure 
plot(x_plot,w)
hold on
plot(x_plot,wl)
legend('nlin','lin')


figure 
plot(x_plot,v)
hold on
plot(x_plot,vl)
legend('nlin','lin')

figure 
plot(x_plot,thx)
hold on
plot(x_plot,thxl)
legend('nlin','lin')

figure 
plot(x_plot,thy)
hold on
plot(x_plot,thyl)
legend('nlin','lin')






