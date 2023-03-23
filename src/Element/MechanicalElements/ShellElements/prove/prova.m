clearvars
close all
clc

%%

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = .1;     % [m] beam's out-of-plane thickness

alpha_T = 11.7e-6; % thermal expansion coefficient

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'THERMAL_EXPANSION_COEFFICIENT', alpha_T);


myElementConstructor = @()Quad8Shell(thickness, myMaterial,2);

% MESH_____________________________________________________________________
Lx = 3;
Ly = .2;
nx = 30;
ny = 5;


% nx = 110;
% ny = 20;

nx = 20;
ny = 3;


eltype = 'QUAD8';
%'QUAD4'

[nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,eltype);


myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

height_midspan = 0.1;

w = height_midspan/Lx;
R = (Lx/4 + Lx*w^2)/(2*w);
yy = R - Lx*w;
xx = Lx/2;
th0 = atan(yy/xx);

x = nodes(:,1);
zNodes = -R*sin(th0)+sqrt(-(x-Lx/2).^2+R^2);

nNodes = myMesh.nNodes;

%set height z
for ii = 1:myMesh.nElements
   
   myMesh.Elements(ii).Object.setZ(zNodes);
    
end

%%
nel = 12;
objp = myMesh.Elements(nel).Object;
X = zeros(nNodes*5,1);
T = 30*randn(nNodes,1);
gradT = 50*randn(nNodes,1);
[T1,Fth] = objp.tangent_stiffness_and_force(X,T,gradT);
M = objp.mass_matrix();

T2 = objp.T2();
T3 = objp.T3();

rng('default');
rng(13);

xtest = 1e-2*randn(nNodes*5,1);
[xtest_el,~,~] = objp.extract_element_data(xtest,T,gradT);
F0 = Fth;
F1 = T1*xtest_el;
F2 = double(ttv(T2,{xtest_el,xtest_el},[2,3]));
F3 = double(ttv(T3,{xtest_el,xtest_el,xtest_el},[2,3,4]));

Ftens = F0 + F1 + F2 + F3;

[~,Fexact] = objp.tangent_stiffness_and_force(xtest,T,gradT);

diff = Fexact-Ftens;
diffcomp = [Fexact-Ftens,Fexact];
compPerc = (Fexact-Ftens)./Fexact*100;
FF = [F0,F1,F2,F3,Fexact,Fexact-Ftens];

%% stiffness derivative (prova)
rng(13);
xtest = 1e-3*randn(nNodes*5,1);
rng(18);
vtest = 1e-2*randn(nNodes*5,1);
Kd = objp.stiffness_derivative(xtest,vtest);
Kinitial = objp.tangent_stiffness_and_force(xtest,T,gradT);

eps = 1e-7;
xperturbed = xtest + eps*vtest;
[Kper_exact,~] = objp.tangent_stiffness_and_force(xperturbed,T,gradT);
dKexact = Kper_exact - Kinitial;
dKper = Kd*eps;

diff_Kd = (dKexact - dKper)./dKexact*100;
max(max(diff_Kd))
% FF3 = objp.F3(xtest,T,gradT);
% FF2 = objp.F2(xtest,T,gradT);
% 
% FFF = Fexact - (F0 + F1 + FF2 + FF3);
% 
% diffF3 = FF3-F3;
% diffF2 = FF2-F2;


%% Boundary conditions
myMesh.set_essential_boundary_condition([nset{1}],1:5,0)
myMesh.set_essential_boundary_condition([nset{3}],1:5,0)

%% Assembly
%T field
T = 0*ones(nNodes,1);
gradT = 0*randn(nNodes,1);

Assembly = Assembly(myMesh);
M = Assembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( myMesh.nDOFs, 1);
[K,~] = Assembly.tangent_stiffness_and_force(u0,T,gradT);

% store matrices
Assembly.DATA.K = K;
Assembly.DATA.M = M;

PlotMesh([nodes,zNodes],elements);

%% Evaluate tensorial forces

% nDOFs = Assembly.Mesh.nDOFs;
% %test case. u0, T and grad T are randomly generated
% u0 = randi(5,nDOFs,1);
% T = randi(25,myMesh.nNodes,1);
% gradT = randi(15,myMesh.nNodes,1);
% 
% tstart = tic;
% 
% %compute constant force and T1
% [K0,Fconst] = Assembly.tangent_stiffness_and_force(0*u0,T,gradT);
% 
% %compute T2
% T2 = Assembly.tensor('T2',[nDOFs, nDOFs, nDOFs], [2,3]);
% F2 = double(ttv(T2,{u0,u0},[2,3]));
% 
% %compute T3
% T3 = Assembly.tensor('T3',[nDOFs, nDOFs, nDOFs, nDOFs], [2,3,4]);
% F3 = double(ttv(T3,{u0,u0,u0},[2,3,4]));
% 
% tend = toc(tstart);
% 
% %approximated force value
% FwithTens = K0*u0 + F2 + F3 + Fconst;
% 
% %compute exact force
% [~,Fexact] = Assembly.tangent_stiffness_and_force(u0,T,gradT);
% diffF = Fexact - FwithTens;
% diffFperc = diffF./Fexact*100;
% 
% 

%% Vibration modes
n_VMs = 2; % first n_VMs modes with lowest frequency calculated 
Kc = Assembly.constrain_matrix(K);
Mc = Assembly.constrain_matrix(M);
[V0,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
V0 = V0(:,ind);
for ii = 1:n_VMs
    V0(:,ii) = V0(:,ii)/max(sqrt(sum(V0(:,ii).^2,2)));
end
V0 = Assembly.unconstrain_vector(V0);

V1 = reshape(V0(:,1),5,[]).';
V2 = reshape(V0(:,2),5,[]).';

figure
PlotFieldonDeformedMesh([nodes,zNodes],elements, [V1(:,1), V1(:,2), V1(:,3)], 'factor',max(nodes(:,2)));
title(['VM1 f: ', num2str(f0(1)),' Hz']);
figure
PlotFieldonDeformedMesh([nodes,zNodes],elements, [V2(:,1), V2(:,2), V2(:,3)], 'factor',max(nodes(:,2)))
title(['VM2 f: ', num2str(f0(2)),' Hz']);

%% Static analysis 
% Nodal force
F = zeros(myMesh.nDOFs,1);
nf = find_node(Lx/2,Ly/2,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
%F(node_force_dofs(2)) = +10e5;
F(node_force_dofs(3)) = +10e5;


%Temperature field
T = 10*ones(nNodes,1);
gradT = 0*randn(nNodes,1);

[K,F0] = Assembly.tangent_stiffness_and_force(u0,T,gradT);

% store matrices
Assembly.DATA.K = K;
Assembly.DATA.F0 = F0;

u_lin = Assembly.solve_system(K, F, F0);
ULIN = reshape(u_lin,5,[]).';	% Linear response
 [~,u] = static_eq(Assembly, F,'vararginTanStiffForce',{T,gradT},'display', 'iter-detailed');
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
node2plot = find_node(plot_dof_location*Lx/3,Ly/2,[],nodes); % node to plot
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






