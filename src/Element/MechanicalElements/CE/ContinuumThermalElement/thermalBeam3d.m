% EXAMPLE: beam meshed with 3D elements
clear;
close all;
clc

% ***** Elements with linear shape functions (very slow mesh convergence, 
%     good for fast code-testing)
%elementType = 'HEX8';
% elementType = 'TET4';

% ***** Elements with quadratic shape functions
% elementType = 'TET10';
%elementType = 'WED15';
elementType = 'HEX20';

%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio

alpha = 11.7e-6;    %thermal expansion coefficient

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu, 'THERMAL_EXPANSION_COEFFICIENT',alpha);
% Element
switch elementType
    case 'HEX8'
        myElementConstructor = @()Hex8Element(myMaterial);
    case 'HEX20'
        myElementConstructor = @()Hex20TElement(myMaterial);
    case 'TET4'
        myElementConstructor = @()Tet4Element(myMaterial);
    case 'TET10'
        myElementConstructor = @()Tet10Element(myMaterial);
    case 'WED15'
        myElementConstructor = @()Wed15Element(myMaterial);
%         quadrature = struct('lin', 5, 'tri', 12);
%         myElementConstructor = @()Wed15Element(myMaterial, quadrature);
end

% MESH_____________________________________________________________________

Lx = 0.4;
Ly = .25;
t = .001;
height_midspan = 0.008;


nx = 30;
ny = 20;
nz = 2;

[nodes, elements, nset] = mesh_3Dparallelepiped(elementType,Lx,Ly,t,nx,ny,nz);
% % Alternatively, one can write an input file in ABAQUS and read it as:
% filename = 'Job-BeamHex';
% [nodes, elements, nset, elset] = mesh_ABAQUSread(filename); % HEX20 mesh


%add to the nodes the vertical displacements due to curvature
w = height_midspan/Lx;
R = (Lx/4 + Lx*w^2)/(2*w);
yy = R - Lx*w;
xx = Lx/2;
th0 = atan(yy/xx);

x = nodes(:,1);
zNodes = -R*sin(th0)+sqrt(-(x-Lx/2).^2+R^2);

nodes(:,3) = nodes(:,3) + zNodes;

%create mesh
myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);
PlotMesh(nodes,elements,0);

% MESH > BOUNDARY CONDITONS
myMesh.set_essential_boundary_condition([nset{1} nset{4}],1:3,0)
myMesh.set_essential_boundary_condition([nset{2} nset{5}],1:3,0)

%myMesh.set_essential_boundary_condition([nset{1}],1:3,0)

% ASSEMBLY ________________________________________________________________
BeamAssembly = Assembly(myMesh);
M = BeamAssembly.mass_matrix();
nNodes = size(nodes,1);

T = 0*ones(myMesh.nNodes,1);
u0 = zeros( myMesh.nDOFs, 1);
[K_cold,~] = BeamAssembly.tangent_stiffness_and_force(u0,T);

% store matrices
BeamAssembly.DATA.K = K_cold;
BeamAssembly.DATA.M = M;


%% vibration modes of cold structure                                      

% Eigenvalue problem_______________________________________________________
n_VMs = 3; % first n_VMs modes with lowest frequency calculated


Kc = BeamAssembly.constrain_matrix(K_cold);
Mc = BeamAssembly.constrain_matrix(M);
[V0,om] = eigs(Kc,Mc,n_VMs,'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
V0 = V0(:,ind);
for ii = 1:n_VMs
    V0(:,ii) = V0(:,ii)/max(sqrt(sum(V0(:,ii).^2,2)));
end
V0 = BeamAssembly.unconstrain_vector(V0);

% PLOT
for mod = [1,2,3]
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes,elements,0);
v1 = reshape(V0(:,mod),3,[]).';
PlotFieldonDeformedMesh(nodes,elements,v1);
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])
end
%drawnow

%% static analysis
F = zeros(myMesh.nDOFs,1);
nf = find_node(Lx/2,Ly/2,t,nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F(node_force_dofs(3)) = +50;

%Temperature field
T = 10*ones(nNodes,1);

%extract tangent stiffness for linear solution
u0 = zeros( myMesh.nDOFs, 1);
[K,gth] = BeamAssembly.tangent_stiffness_and_force(u0,T);

u_lin = BeamAssembly.solve_system(K, F-gth);
ULIN = reshape(u_lin,3,[]).';	% Linear response
[~,u] = static_equilibrium_thermal(BeamAssembly, F, T, 'display', 'iter-detailed');
UNL = reshape(u,3,[]).';        % Nonlinear response

%linear solution
Ulinx = ULIN(:,1);
Uliny = ULIN(:,2);
Ulinz = ULIN(:,3);

%nonlinear solution
UNlinx = UNL(:,1);
UNliny = UNL(:,2);
UNlinz = UNL(:,3);


figure
PlotFieldonDeformedMesh(nodes,elements, [Ulinx, Uliny, Ulinz], 'factor',20)
hold on;
PlotMesh(nodes,elements);
grid on
title('linear static solution');
axis

figure
PlotFieldonDeformedMesh(nodes,elements, [UNlinx, UNliny, UNlinz], 'factor',20)
hold on;
PlotMesh(nodes,elements);
title('non linear static solution');
grid on
axis

%% analyze equilibrium point

[Keq,Feq] = BeamAssembly.tangent_stiffness_and_force(u,T);

n_VMs = 3; % first n_VMs modes with lowest frequency calculated 
Keq_c = BeamAssembly.constrain_matrix(Keq);
[V0,om_eq] = eigs(Keq_c, Mc, n_VMs, 'SM');
[f0_eq,ind] = sort(sqrt(diag(om_eq))/2/pi);

V0 = V0(:,ind);
for ii = 1:n_VMs
    V0(:,ii) = V0(:,ii)/max(sqrt(sum(V0(:,ii).^2,2)));
end
V0 = BeamAssembly.unconstrain_vector(V0);

VMs = cell(n_VMs,1);
for ii = 1:n_VMs  
    VMs{ii,1} = reshape(V0(:,ii),3,[]).';    
end

%plot VMs
for ii = 1:n_VMs
    figure; hold on;
    VM2plot = VMs{ii};
    PlotFieldonDeformedMesh(nodes,elements, [VM2plot(:,1),...
        VM2plot(:,2), VM2plot(:,3)], 'factor',max(nodes(:,2)));
    title(['VM at eq ',num2str(ii),', f: ', num2str(f0_eq(ii)),' Hz']);
end

du = 1e-9*randn(length(u),1);
[~,F1] = Assembly.tangent_stiffness_and_force(u+du,T,gradT);

F1app = Feq + Keq*du;

diffF = norm(F1-F1app)/norm(F1-Feq)*100;

diffF;

dd = (F1-F1app)./(F1-Feq)*100;

ddmax = max(dd);


