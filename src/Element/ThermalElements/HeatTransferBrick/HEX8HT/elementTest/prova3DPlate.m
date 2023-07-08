close all
clc
clearvars

%% test of thermal fem

% DATA ___________________________________________________________________
E = 200e9;              % 70e9 % 200e9 % Young's modulus
rho = 7850;             % 2700 % 7850 % density
k = 45;                  % conductivity
c = 420;                % specific heat

myMaterial  = Material();
set(myMaterial,'THERMAL_CONDUCTIVITY', k, 'SPECIFIC_HEAT', c, 'DENSITY', rho); 


% Element
myElementConstructor = @()HEX8ElementHT(myMaterial);

% MESH_____________________________________________________________________
l = 2;
w = 2;
t = .9;
nx = 15;
ny = 15;
nz = 15;
elementType = 'HEX8';
[nodes, elements, nset]=mesh_3Dparallelepiped_mod(elementType,l,w,t,nx,ny,nz);

%nodes(:,3) = nodes(:,3)+1*t*(cos(2*pi/l*nodes(:,1))+1);
%PlotMesh(nodes,elements,1);

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

%%
% BOUNDARY CONDITONS
Timposed = 0;
Timposed1 = 0;%250;
Timposed2 = 0;%150;
setNatBCs = [nset.nodes{1}];
%setNatBCs = [nset.nodes{1},nset.nodes{2},nset.nodes{5}];
setNatBCs1 = nset.nodes{4};
setNatBCs2 = nset.nodes{5};

myMesh.set_essential_boundary_condition(setNatBCs,1,Timposed)
myMesh.set_essential_boundary_condition(setNatBCs1,1,Timposed1)
myMesh.set_essential_boundary_condition(setNatBCs2,1,Timposed2)

% NATURAL BOUNDARY CONDITIONS (CREATE NEUMANN ELEMENTS)
setNodeNeu_index = 6;
nodesNeumann = nset.nodes{setNodeNeu_index};
nodesCoordNeumann = nodes(nset.nodes{setNodeNeu_index},:);
elementsNeumann =  nset.connectivity{setNodeNeu_index};
myElementConstructor = @()HEX8NeuElementHT(myMaterial);
myMesh.create_elements_table(elementsNeumann,myElementConstructor,'isBoundary',true);

% ASSEMBLY OF SYSTEM MATRICES
Assembly = ThermalAssembly(myMesh);

K = Assembly.conductivity_matrix();
M = Assembly.thermal_capacity_matrix();
LM = Assembly.load_matrix();
QB = Assembly.thermal_flux_matrix();

TfromBCs = zeros(Assembly.Mesh.nNodes,1);
TfromBCs(setNatBCs) = Timposed;
TfromBCs(setNatBCs1) = Timposed1;
TfromBCs(setNatBCs2) = Timposed2;
loadNaturalBCs = Assembly.constrain_vector(-K*TfromBCs); %thermal load from imposed BCs

%external applied thermal load (flux)
q = zeros(Assembly.Mesh.nNodes,1);
%q(nodesNeumann) = -200;



%externally applied thermal load (on shell surface)
xc = l/2;
yc = w/2;
d = norm([xc,yc])^2;
WattsInSurface = 400;
xnodes = nodes(:,1);
ynodes = nodes(:,2);
coeffW = WattsInSurface/l/w;
qS = (-((xnodes-xc).^2+(ynodes-yc).^2)+d)*1/d*coeffW;
q(nodesNeumann) = qS(nodesNeumann);




%constrain thermal conductivity matrix and solve system
K_c = Assembly.constrain_matrix(K);
qLoad = Assembly.constrain_vector(QB*q);

Tsol = K_c\(qLoad+loadNaturalBCs);

Tsol = Assembly.unconstrain_vector(Tsol);


%PLOT OF RESULTS
figure
PlotFieldonMesh(nodes,elements,Tsol)


dl = l/nx;
ypath = w/2;
zpath = t;
nodesPath = zeros(nx+1,1);
for ii = 1:nx+1
    lii = (ii-1)*dl;
    nodesPath(ii) = find_node(lii,ypath,zpath,nodes);
end
xnodesPath = nodes(nodesPath,1);
Tpath = Tsol(nodesPath);

figure
plot(xnodesPath,Tpath,'-x');
title('solution on path along x')
xlabel('x [m]')
ylabel('T [K]')


    

