clearvars
close all
clc

%%

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = .1;     % [m] beam's out-of-plane thickness
alpha_T = 11.7e-7; % thermal expansion coefficient

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'THERMAL_EXPANSION_COEFFICIENT', alpha_T);


myElementConstructor = @()Quad8Shell(thickness, myMaterial);

% MESH_____________________________________________________________________
Lx = 1;
Ly = 1;
nx = 1;
ny = 1;

[nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,'QUAD8');


myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

%%
objp = myMesh.Elements(1).Object;
X = zeros(620,1);
T = zeros(124,1);
gradT = zeros(124,1);
[K,F] = objp.tangent_stiffness_and_force(X,T,gradT);

%% Boundary conditions
 myMesh.set_essential_boundary_condition([1,4],1:5,0)

%% Assembly
%T field
T = 100*ones(8,1);
gradT = 0*randn(8,1);

Assembly = Assembly(myMesh);
M = Assembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( myMesh.nDOFs, 1);
[K,~] = Assembly.tangent_stiffness_and_force(u0,T,gradT);
K = full(K);
M = full(M);

isKsymm = K-K';
isMsymm = M-M';

% store matrices
Assembly.DATA.K = K;
Assembly.DATA.M = M;

PlotMesh(nodes,elements);

Kc = Assembly.constrain_matrix(K);
Mc = Assembly.constrain_matrix(M);

%% Static analysis 
% Nodal force
F = zeros(myMesh.nDOFs,1);
nf = 2; % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );

F(node_force_dofs(5)) = +10e5;
F(node_force_dofs(5)) = 0;

u_lin = Assembly.solve_system(K, F);
u_lin = reshape(u_lin,5,[]);


% figure 
% hold on
% scl = 20;
% for ii = 1:4
% plot(nodes(ii,1),nodes(ii,2),'x','color','b','markersize',10);
% plot(nodes(ii,1)+scl*u_lin(1,ii),nodes(ii,2)+scl*u_lin(2,ii),'x','color','r','markersize',10)
% axis([-0.1,1.1,-0.1,1.1])
% end
% grid on
% axis([-0.1,1.1,-0.1,1.1])
% 
% w = u_lin(3,:);
% thx = u_lin(4,:);
% thy = u_lin(5,:);

%nonlinear
[~,unl] = static_eq_shell_thermal(Assembly, F, T, gradT, 'display', 'iter-detailed');
u_nl = reshape(unl,5,[]).'; 

figure 
hold on
scl = 20;
for ii = 1:4
plot(nodes(ii,1),nodes(ii,2),'x','color','b','markersize',10);
plot(nodes(ii,1)+scl*u_nl(ii,1),nodes(ii,2)+scl*u_nl(ii,2),'x','color','r','markersize',10)
axis([-0.1,1.1,-0.1,1.1])
end
grid on
axis([-0.1,1.1,-0.1,1.1])


%tangent stiffness
[Ktg,Feq] = Assembly.tangent_stiffness_and_force(unl,T,gradT);

du = 1e-3*randn(40,1);

u_new = unl+du;

dfapp = Ktg*du;

[~,Fnew] = Assembly.tangent_stiffness_and_force(u_new,T,gradT);

err = 100*(Feq+dfapp - Fnew)./(Fnew-Feq);






return

ULIN = reshape(u_lin,5,[]).';	% Linear response
u = static_eq_shell_thermal(Assembly, F, T, gradT, 'display', 'iter-detailed');
UNL = reshape(u,5,[]).';        % Nonlinear response


% 
% node_plot_dofs = get_index(nset{2}, myMesh.nDOFPerNode );
% dofs_u = node_plot_dofs(1:5:end);
% dofs_v = node_plot_dofs(2:5:end);
% dofs_w = node_plot_dofs(3:5:end);
% dofs_thx = node_plot_dofs(4:5:end);
% dofs_thy = node_plot_dofs(5:5:end);

x_plot = (nodes(nset{2}))';

u = UNL(nset{2},1);
v = UNL(nset{2},2);
w = UNL(nset{2},3);
thx = UNL(nset{2},4);
thy = UNL(nset{2},5);

ul = ULIN(nset{2},1);
vl = ULIN(nset{2},2);
wl = ULIN(nset{2},3);
thxl = ULIN(nset{2},4);
thyl = ULIN(nset{2},5);


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






