% author: ALEXANDER SACCANI - PHD CANDIDATE, ETH ZURICH

% verify natural frequencies of hot temperature beam
% the model considered is a clamped-clamped beam with different
% temperatures distributions:
% 1) constant T
% 2) ramp
% 3) parabola

close all
clc
clearvars

model = 'arc'; %valid options are : 'straightbeam' or 'arc'

if strcmp(model,'straightbeam')
    
    %% straight beam model construction

    %GEOMETRY_____________________________________________________________
    l = 0.1; 
    h = 1e-3;
    b = 1e-2; 

    %MATERIAL______________________________________________________________
    E       = 70e9;  % 70e9 % 200e9 % Young's modulus
    rho     = 2700; % 2700 % 7850 % density
    nu      = 0.3;    % nu
    kappa   = 1e8; % material damping modulus 1e8
    alpha_T = 23.1e-6; % thermal expansion coefficient

    myThermalBeamMaterial = KirchoffMaterial();
    set(myThermalBeamMaterial,'YOUNGS_MODULUS', E,'DENSITY',rho,...
        'POISSONS_RATIO',nu,'THERMAL_EXPANSION_COEFFICIENT', alpha_T);

    %MESH__________________________________________________________________
    nElements = 60;

    %define nodes
    dx = l/nElements; %x spacing between nodes

    nodes_x = (0:dx:l).'; %x coordinates for nodes
    nodes_y = zeros(length(nodes_x),1); %y coordinates for nodes

    nNodes = length(nodes_x); %number of nodes
    nDofPerNode = 3;

    %create mesh from nodes
    Nodes = [nodes_x, nodes_y];
    BeamMesh = Mesh(Nodes); 

    %define elements 
    Elements = [1:nNodes-1;2:nNodes].'; % ROW -> # element, 
                                        % COLUMN -> # of node connected
    myElementConstructor = @()ThermalBeamElement(b, h, myThermalBeamMaterial); %type of element

    %create element table from elements
    BeamMesh.create_elements_table(Elements,myElementConstructor);


    %BOUNDARY CONDITIONS_______________________________________________________
    BeamMesh.set_essential_boundary_condition([1 BeamMesh.nNodes],[1 2 3],0) % Doubly clamped beam
    %BeamMesh.set_essential_boundary_condition(1,[1 2 3],0) % clamped beam  (node,[constrained dofs],imposed disp.)

    BeamAssembly = Assembly(BeamMesh);

    nDofsF = BeamAssembly.Mesh.EBC.nDOFs; % #dofs free (without bc)
    nDofsC = nDofsF - size(BeamAssembly.Mesh.EBC.constrainedDOFs,1); % #dofs constrained (with bc)

end

if strcmp(model,'arc')
    
    %% arc model construction

    %GEOMETRY__________________________________________________________________
    l = 0.1; 
    h = 1e-3;
    b = 1e-2; 
    height_midspan = 1e-2;

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

    %plot nodes
    figure; title('beam mesh')
    plot(nodes_x,nodes_y,'*-');
    axis('equal')

    %create mesh from nodes
    Nodes = [nodes_x, nodes_y];
    BeamMesh = Mesh(Nodes); 

    %define elements (connectivity matrix)
    Elements = [1:nNodes-1;2:nNodes].'; % ROW -> # element, 
                                        % COLUMN -> # of node connected
    myElementConstructor = @()ThermalBeamElement(b, h, myThermalBeamMaterial); %type of element

    %create element table from elements
    BeamMesh.create_elements_table(Elements,myElementConstructor);


    %BOUNDARY CONDITIONS_______________________________________________________
    BeamMesh.set_essential_boundary_condition([1 BeamMesh.nNodes],[1 2 3],0) % Doubly clamped beam
    %BeamMesh.set_essential_boundary_condition(1,[1 2 3],0) % clamped beam  (node,[constrained dofs],imposed disp.)

    BeamAssembly = Assembly(BeamMesh);

    nDofsF = BeamAssembly.Mesh.EBC.nDOFs; % #dofs free (without bc)
    nDofsC = nDofsF - size(BeamAssembly.Mesh.EBC.constrainedDOFs,1); % #dofs constrained (with bc)

end

%GET USEFUL QUANTITIES FROM THE CREATED MODEL______________________________
K_cold = BeamAssembly.stiffness_matrix(); % tangent stiffness matrix for T = 0, u = 0
M = BeamAssembly.mass_matrix();
C = BeamAssembly.damping_matrix();

K_C = BeamAssembly.constrain_matrix(K_cold); % should be tangent stiffness matrix for T = 0, u = 0
M_C = BeamAssembly.constrain_matrix(M);

%% Modal Analysis of the cold structure

n_VMs = 3; % first n_VMs modes with lowest frequency calculated 

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


%% Thermal VMs analysis 

% %constant T
% T_dyn_p = @(T_ampl) T_ampl; % field parametrized with thermal amplitude
% p_sampl = [0;5;10;30]; 

% %ramp
% T_dyn_p = @(T_ampl) T_ampl/l*nodes_x; %field parametrized with thermal amplitude at right end
% p_sampl = [5;10;15;20]; 

% %parabola
% T_dyn_p = @(T_ampl) -4*T_ampl/l^2*nodes_x.^2 + 4*T_ampl/l*nodes_x; %field parametrized with thermal amplitude at midspan
% p_sampl = [0;20]; 

% % sine pulse
% T_ampl = 100; %width of thermal pulse
% p = 3e-2; %width of sine pulse
% x0 = @(xc) xc - p/2; %left extreme of pulse
% T_dyn_p = @(xc) T_ampl*(sin(pi*(nodes_x-x0(xc))/p)).^2.*(heaviside(nodes_x-x0(xc)) - ...
%     heaviside((nodes_x-x0(xc))-p)); %define the temperature as a function of xc only (T profile parametrized with xc
% %p_sampl = [-p/2;0;0.01;0.05;0.085;0.1]; % field parametrized with xc (centre of pulse)
% nsampls = 4;
% p_sampl = [-p/2,l/nsampls:l/nsampls:l].';

%half beam cold - half hot
perc_hot = 80;
T_dyn_p = @(T_ampl) T_ampl*(1 - heaviside(nodes_x - (perc_hot/100)*l));
p_sampl = [0,5,10,50];




T_sampl = zeros(nNodes,length(p_sampl)); % T nodal distribution samples

%generate corresponding temperature profiles
for ii = 1:length(p_sampl)
    T_sampl(:,ii) = T_dyn_p(p_sampl(ii));
end

VMs_at_eq = 1; %do you want to compute VMs at thermal equilibrium? (set to 1 if yes, otherwise set to 0)
[VMs,static_eq,omega] = VM_thermal(BeamAssembly,T_sampl,VMs_at_eq,5);

%plot VMs__________________________________________________________________
color_list = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE'};
T_sampl_2_plot = 1:length(p_sampl);
mode2plot = [1,2,3,4,5];

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
        plot(nodes_x + u_ii, nodes_y + v_ii,'.-','markersize',10,'color',color_list{jj}); hold on;
        plot(nodes_x - u_ii, nodes_y - v_ii,'.-','markersize',10,'color',color_list{jj}); hold on;
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

disp('natural frequencies are [cycles/s]: (row -> VM, column -> T distribution)');
disp(omega/(2*pi));