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
height_midspan = 4e-3; %height of midspan 

%Mesh______________________________________________________________________
nx = 15; % # elements along x
ny = 8; % # elements along y

eltype = 'QUAD8'; %'QUAD4'
myElementConstructor = @()ThermQuad8Shell(thickness, myMaterial,2);
[nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,eltype);

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

nNodes = myMesh.nNodes;

nDofPerNode = 5;

elementsPlot = elements(:,1:4); %for plots

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
figure
PlotMesh([nodes,zNodes],elementsPlot);
axis on

%Boundary conditions: all sides are clamped________________________________
myMesh.set_essential_boundary_condition([nset{1}],1:3,0)
myMesh.set_essential_boundary_condition([nset{2}],1:3,0)
myMesh.set_essential_boundary_condition([nset{3}],1:3,0)
myMesh.set_essential_boundary_condition([nset{4}],1:3,0)

%Assembly__________________________________________________________________
Assembly = Assembly(myMesh);
M = Assembly.mass_matrix();
u0 = zeros( myMesh.nDOFs, 1);

%Compute and store Mass matrix, K cold_____________________________________
%T field (zeros)
T = zeros(nNodes,1);
gradT = zeros(nNodes,1);
[Kc,~] = Assembly.tangent_stiffness_and_force(u0,T,gradT);

alpha = 0.01;
beta = 0.01;
D = alpha*Kc + beta*M;

% store matrices
Assembly.DATA.Kc = Kc; %Kcold
Assembly.DATA.M = M; %mass matrix
Assembly.DATA.C = D; %damping matrix


%% Vibration modes (cold)
% nVMs = 5; % first n_VMs modes with lowest frequency calculated 
% 
% Kconstr = Assembly.constrain_matrix(Assembly.DATA.Kc);
% Mconstr = Assembly.constrain_matrix(Assembly.DATA.M);
% 
% [VMsCold,om] = eigs(Kconstr, Mconstr, nVMs, 'SM');
% [fcold,ind] = sort(sqrt(diag(om))/2/pi);
% 
% VMsCold = VMsCold(:,ind);
% for ii = 1:nVMs
%     VMsCold(:,ii) = VMsCold(:,ii)/max(sqrt(sum(VMsCold(:,ii).^2,2)));
% end
% VMsCold = Assembly.unconstrain_vector(VMsCold);
% 
% %plot VMs__________________
% VMs = cell(nVMs,1);
% for ii = 1:nVMs  
%     VMs{ii,1} = reshape(VMsCold(:,ii),5,[]).';    
% end
% 
% for ii = 1:nVMs
%     figure; hold on;
%     VM2plot = VMs{ii};
%     PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [VM2plot(:,1),...
%         VM2plot(:,2), VM2plot(:,3)], 'factor',max(nodes(:,2)));
%     title(['VM ',num2str(ii),', f: ', num2str(fcold(ii)),' Hz']);
% end

%% VMs and MDs Cold

T_ampls = 0;

n_sampl = length(T_ampls);
T_samples = zeros(myMesh.nNodes,n_sampl);
gradT_samples = zeros(myMesh.nNodes,n_sampl);

for ii = 1:n_sampl
    T_samples(:,ii) = ones(myMesh.nNodes,1)*T_ampls(ii);
end

nVMs = 5;

logic_MD = 1;
VMMD = VMMD_thermal(Assembly,nVMs,logic_MD,T_samples,gradT_samples);
fcold = VMMD.omega/2/pi;
VMsCold = VMMD.VMs;
MDsCold = VMMD.MDs.data;
MDsNames = VMMD.MDs.names;
eqPoint = VMMD.eq;

%Plot VMs
for jj = 1:length(T_ampls)
    
    VMs = cell(nVMs,1);
    for ii = 1:nVMs  
        VMs{ii,1} = reshape(VMsCold{jj}(:,ii),5,[]).';    
    end

    for ii = 1:nVMs
        figure; hold on;
        VM2plot = VMs{ii};
        PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [VM2plot(:,1),...
            VM2plot(:,2), VM2plot(:,3)], 'factor', 15*max(nodes(:,2)));
        title(['VM ',num2str(ii),' dT = ',num2str(T_ampls(jj)),'K, f: ', num2str(fcold(ii,jj)),' Hz']);
    end

end


%plot MDs
nMDs = size(MDsCold{1,1},2);

for jj = 1:length(T_ampls)
    
    MDs = cell(nMDs,1);
    for ii = 1:nMDs  
       MDs{ii,1} = reshape(MDsCold{jj}(:,ii),5,[]).';    
    end

    for ii = 1:nMDs
        figure; hold on;
        MD2plot = MDs{ii};
        PlotFieldonDeformedMesh([nodes,zNodes],elementsPlot, [MD2plot(:,1),...
            MD2plot(:,2), MD2plot(:,3)], 'factor',1*max(nodes(:,2)));
        title(['MD ',num2str(MDsNames(ii,1)),num2str(MDsNames(ii,2)),'  dT = ',num2str(T_ampls(jj)),'K']);
    end

end
