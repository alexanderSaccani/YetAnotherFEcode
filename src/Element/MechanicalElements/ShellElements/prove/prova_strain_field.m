clearvars
close all
clc

%%

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = .01;     % [m] beam's out-of-plane thickness
alpha_T = 11.7e-7; % thermal expansion coefficient

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'THERMAL_EXPANSION_COEFFICIENT', alpha_T);


myElementConstructor = @()Quad8Shell(thickness, myMaterial);

% MESH_____________________________________________________________________
Lx = 1;
Ly = 0.3;
nx = 30;
ny = 8;

[nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,'QUAD8');


myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

elements2plot = elements(:,1:4);
PlotMesh(nodes,elements2plot,0) 

%% Boundary conditions
myMesh.set_essential_boundary_condition(nset{1},1:5,0);
% myMesh.set_essential_boundary_condition(nset{3},1:5,0);

%% Assembly
Assembly = Assembly(myMesh);

%tang stiff 
M = Assembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( myMesh.nDOFs, 1);
[K,~] = Assembly.tangent_stiffness_and_force(u0,zeros(myMesh.nNodes,1),zeros(myMesh.nNodes,1));

% store matrices
Assembly.DATA.K = K;
Assembly.DATA.M = M;

Kc = Assembly.constrain_matrix(K);
Mc = Assembly.constrain_matrix(M);

%% Static analysis 
% Nodal force
F0 = zeros(myMesh.nDOFs,1);
nf = find_node(Lx,Ly/2,[],myMesh.nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F0(node_force_dofs(2)) = +150;

% %distr load
% node_force_dofs = get_index(nset{3}, myMesh.nDOFPerNode );
% node_force_dofs = reshape(node_force_dofs,5,[]);
% node_force_dofs = node_force_dofs(1,:); %load applied along x
% F0(node_force_dofs') = 120;

T = 0*ones(myMesh.nNodes,1);
gradT = zeros(myMesh.nNodes,1);

[~,uStaticF] = static_eq_shell_thermal(Assembly, F0, T,gradT, 'display', 'iter-detailed');


method = 'membrane_forces';
% nodeSet = 'all';
nodeSet = 'all';
nOut = 3;
field = Assembly.get_field(elements,nodeSet,method,nOut,uStaticF,T,gradT);

%%
elements2plot = elements(:,1:4);

figure
PlotFieldonMesh(nodes,elements2plot ,field(:,1));
title('Nx')

figure
PlotFieldonMesh(nodes,elements2plot ,field(:,2));
title('Ny')

figure
PlotFieldonMesh(nodes,elements2plot ,field(:,3));
title('Nxy')

disp = reshape(uStaticF,5,[]).';
disp = disp(:,1:3);
figure
PlotFieldonDeformedMesh([nodes,zeros(size(nodes,1),1)],elements2plot,disp);
title('disp')

nodeset1 = zeros(ny+1,1);
for ii = 1:ny+1
nf = find_node(Lx/nx*18,Ly/ny*(ii-1),[],myMesh.nodes); % node where to put the force
nodeset1(ii) = nf;
end
xnodeset1 = nodes(nodeset1,1);
ynodeset1 = nodes(nodeset1,2);

nset2plot = nodeset1;

% indNodeSetLx = get_index(nset2plot, myMesh.nDOFPerNode );
% indNodeSetLx = reshape(indNodeSetLx,5,[]);
% indNodeSetLxX = indNodeSetLx(1,:).';
% indNodeSetLxY = indNodeSetLx(2,:).';
%xNodeSet = nodes(nset{3},1);
yNodeSetLx = nodes(nset2plot,2);

figure
plot(yNodeSetLx,disp(nset2plot,1),'*');
title('disp along x')

method = 'membrane_forces';
% nodeSet = 'all';
nOut = 3;
fieldLx = Assembly.get_field(elements,nset2plot,method,nOut,uStaticF,T,gradT);

figure 
plot(yNodeSetLx,fieldLx(:,1),'x')
hold on
plot(yNodeSetLx,field(nset2plot,1),'--')
title('membrane Field X')

figure 
plot(yNodeSetLx,fieldLx(:,2),'x')
hold on
plot(yNodeSetLx,field(nset2plot,2),'--')
title('membrane Field Y')

figure 
plot(yNodeSetLx,fieldLx(:,3),'x')
hold on
plot(yNodeSetLx,field(nset2plot,3),'--')
title('membrane Field XY')

return
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






