close all
clc
clearvars

%% test of thermal fem

% DATA ___________________________________________________________________
E = 200e9;              %70e9 % 200e9 % Young's modulus
rho = 7850;             %2700 % 7850 % density
k = 45;                 %conductivity
c = 420;                %specific heat
convCoeff = 10;         %convection coefficient (numero a caso)
Tsink = 180;            %temperature of the sink
convection = true;      %convection on or off

myMaterial  = Material();
set(myMaterial,'THERMAL_CONDUCTIVITY', k, 'SPECIFIC_HEAT', c, 'DENSITY', rho); 


% MESH_____________________________________________________________________
Lx = 2;
Ly = 1;
thickness = .001;
height_midspan = 1e-10*0.5;
nx = 20;
ny = 20;
elementType = 'QUAD8';
[nodes, elements, bset, belements ] =  mesh_2Drectangle_mod(Lx,Ly,nx,ny,elementType);


%define z to assign to nodes
w = height_midspan/Lx;
R = (Lx/4 + Lx*w^2)/(2*w);
yy = R - Lx*w;
xx = Lx/2;
th0 = atan(yy/xx);
x = nodes(:,1);
zNodes = -R*sin(th0)+sqrt(-(x-Lx/2).^2+R^2);

nodes = [nodes,zNodes];


% Element
myElementConstructor = @()QUAD8ShellHT(myMaterial,thickness);

%nodes(:,3) = nodes(:,3)+1*t*(cos(2*pi/l*nodes(:,1))+1);
%PlotMesh(nodes,elements,1);

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

%%
% BOUNDARY CONDITONS
Timposed = 50;
Timposed1 = 0;
Timposed2 = 0;%150;
setNatBCs = [bset{1}];
%setNatBCs = [bset.nodes{1},bset.nodes{2},bset.nodes{5}];
setNatBCs1 = bset{3};
setNatBCs2 = bset{4};

myMesh.set_essential_boundary_condition(setNatBCs,1,Timposed)
myMesh.set_essential_boundary_condition(setNatBCs1,1,Timposed1)
myMesh.set_essential_boundary_condition(setNatBCs2,1,Timposed2)

% NATURAL BOUNDARY CONDITIONS (CREATE NEUMANN ELEMENTS)
setNodeNeu_index = 2;
nodesNeumann = bset{setNodeNeu_index};
nodesCoordNeumann = nodes(bset{setNodeNeu_index},:);
elementsNeumann =  belements{setNodeNeu_index};
myElementConstructor = @()QUAD8ShellNeuHT(myMaterial,thickness);
myMesh.create_elements_table(elementsNeumann,myElementConstructor,'isBoundary',true);

% ASSEMBLY OF SYSTEM MATRICES
Assembly = ThermalAssembly(myMesh);

K = Assembly.conductivity_matrix();
M = Assembly.thermal_capacity_matrix();
QS = Assembly.load_matrix();
QB = Assembly.thermal_flux_matrix();

TfromBCs = zeros(Assembly.Mesh.nNodes,1);
TfromBCs(setNatBCs) = Timposed;
TfromBCs(setNatBCs1) = Timposed1;
TfromBCs(setNatBCs2) = Timposed2;
loadNaturalBCs = Assembly.constrain_vector(-K*TfromBCs); %thermal load from imposed BCs

%externally applied thermal load (flux, on boundary)
qB = zeros(Assembly.Mesh.nNodes,1);
WattsIn = -60;%-60;
qPerUnitSurface = WattsIn/thickness/Lx;
qB(nodesNeumann) = qPerUnitSurface;

%externally applied thermal load (on shell surface)
xc = Lx/2;
yc = Ly/2;
d = norm([xc,yc])^2;
WattsInSurface = 10;
xnodes = nodes(:,1);
ynodes = nodes(:,2);
coeffW = WattsInSurface/Lx/Ly;
qS = (-((xnodes-xc).^2+(ynodes-yc).^2)+d)/d*coeffW;

%qS = zeros(Assembly.Mesh.nNodes,1) +20;

%constrain thermal conductivity matrix and solve system
K_c = Assembly.constrain_matrix(K);
a = min(QB*qB);
b = min(QS*qS);
qLoad = Assembly.constrain_vector(QB*qB+QS*qS);

Tsol = K_c\(qLoad+loadNaturalBCs);

Tsol = Assembly.unconstrain_vector(Tsol);


%PLOT OF RESULTS
elementsPlot = elements(:,1:4); %for plots
figure
PlotFieldonMesh(nodes,elementsPlot,Tsol)

%path along x
dl = Lx/nx;
ypath = 0;
nodesPathx = zeros(nx+1,1);
for ii = 1:nx+1
    lii = (ii-1)*dl;
    nodesPathx(ii) = find_node(lii,ypath,0,[nodes(:,1:2),zeros(size(nodes,1),1)]);
end
xnodesPath = nodes(nodesPathx,1);
Tpathx = Tsol(nodesPathx);

figure
plot(xnodesPath,Tpathx,'-x');
title('solution on path along x')
xlabel('x [m]')
ylabel('T [K]')


%path along y
dl = Ly/nx;
xpath = Lx/2;
nodesPathy = zeros(ny+1,1);
for ii = 1:ny+1
    lii = (ii-1)*dl;
    nodesPathy(ii) = find_node(xpath,lii,0,[nodes(:,1:2),zeros(size(nodes,1),1)]);
end
ynodesPath = nodes(nodesPathy,2);
Tpathy = Tsol(nodesPathy);

figure
plot(ynodesPath,Tpathy,'-x');
title('solution on path along y')
xlabel('x [m]')
ylabel('T [K]')


    

