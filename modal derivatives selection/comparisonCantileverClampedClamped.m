close all
clearvars
clc
%% Construct model

%Material _________________________________________________________________
E       = 71e9;     % Young's modulus [Pa]
rho     = 2795;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 

myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);

%Geometry__________________________________________________________________
% Lx = 0.65;
% Ly = .125;
% thickness = .003;     % thickness of plate
% height_midspan = 5e-10; %height of midspan 

Lx = 0.55;
Ly = .06;
thickness = .00225;     % thickness of plate
height_midspan = 5e-10; %height of midspan 

%Mesh______________________________________________________________________
nx = 15; % # elements along x
ny = 8; % # elements along y

eltype = 'QUAD8'; %'QUAD4'
myElementConstructor = @()Quad8Shell(thickness, myMaterial,2);
[nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,eltype);

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

nNodes = myMesh.nNodes;

nDofPerNode = 5;

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
elementsPlot = elements(:,1:4); %for plots
PlotMesh([nodes,zNodes],elementsPlot);

%Boundary conditions: all sides are clamped________________________________
myMesh.set_essential_boundary_condition([nset{1}],1:5,0)
%myMesh.set_essential_boundary_condition([nset{3}],1:5,0)

%Assembly__________________________________________________________________
Assembly = Assembly(myMesh);
M = Assembly.mass_matrix();
u0 = zeros( myMesh.nDOFs, 1);

%Compute and store Mass matrix, K matrix
[K,~] = Assembly.tangent_stiffness_and_force(u0);

% store matrices
Assembly.DATA.K = K; %K
Assembly.DATA.M = M; %mass matrix


%% VMs analysis 

nVMs = 5; % first n_VMs modes with lowest frequency calculated 

Kconstr = Assembly.constrain_matrix(Assembly.DATA.K);
Mconstr = Assembly.constrain_matrix(Assembly.DATA.M);

[VMs,om] = eigs(Kconstr, Mconstr, nVMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);

VMs = VMs(:,ind);
for ii = 1:nVMs
    VMs(:,ii) = VMs(:,ii)/max(sqrt(sum(VMs(:,ii).^2,2)));
end
VMs = Assembly.unconstrain_vector(VMs);

%plot VMs 
VMsPlot = cell(nVMs,1); 
for ii = 1:nVMs
    VMsPlot{ii,1} = reshape(VMs(:,ii),5,[]).';
end

VMs2Plot = 3;%1:5; 
for ii = 1:length(VMs2Plot)
    figure; hold on; VM2plotii = VMsPlot{VMs2Plot(ii)};
    PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot,[VM2plotii(:,1),...
        VM2plotii(:,2), VM2plotii(:,3)], 'factor',1.5*max(nodes(:,2)));
    title(['VM ',num2str(VMs2Plot(ii)),', f: ',num2str(f0(VMs2Plot(ii))),' Hz']);
end

%% Rayleigh Damping
r1 = 0.01; %damping of first mode
r2 = 0.01;  %damping of second mode

r = [r1;r2];
om1 = 2*pi*f0(1);
om2 = 2*pi*f0(2);

dampParam = [1 om1^2; 1 om2^2]\r;

C = dampParam(1)*M + dampParam(2)*K;

Assembly.DATA.C = C;


% %% Static Analysis
% % Nodal force
% F = zeros(myMesh.nDOFs,1);
% nf = find_node(0.4*Lx,0.25*Ly,[],nodes); % node where to put the force
% node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
% F(node_force_dofs(3)) = 75;
% 
% 
% [ulin,unl] = static_equilibrium(Assembly, F, 'display', 'iter-detailed');
% UNL = reshape(unl,5,[]).';        % Nonlinear response
% ULIN = reshape(ulin,5,[]).';	% Linear response
% 
% 
% % PLOT
% elementPlot = elements(:,1:4);
% figure('units','normalized','position',[.2 .1 .6 .8])
% scale = 1;
% PlotMesh(nodes, elementPlot, 0);
% PlotFieldonDeformedMesh([nodes,zNodes],elementPlot,[ULIN(:,1),ULIN(:,2),ULIN(:,3)],'factor',scale,'color','k');
% hold on
% PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot,[UNL(:,1),UNL(:,2),UNL(:,3)],'factor',scale,'color','k');
% colormap jet
% title(['LINEAR and NONLINEAR STATIC RESPONSE, scale factor: ' num2str(scale) 'x)']);

%% Modal Derivatives

%Static Modal Derivatives
%how many VMs
nVMs4MDs = 4;
[SMDs, names] = modal_derivatives(Assembly, elements, VMs(:,1:nVMs4MDs));

%Plot static Modal Derivatives
nMDs = size(SMDs,2);
MDs = cell(nMDs,1);
for ii = 1:nMDs  
   MDs{ii,1} = reshape(SMDs(:,ii),5,[]).';    
end

for ii = 1:nMDs
    figure; hold on;
    MD2plot = MDs{ii};
    PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [MD2plot(:,1),...
        MD2plot(:,2), MD2plot(:,3)], 'factor',1*max(nodes(:,2)));
    title(['MD ',num2str(names(ii,1)),num2str(names(ii,2))]);
end


