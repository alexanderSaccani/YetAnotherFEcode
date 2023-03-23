function [VM,static_sol,omega] = VM_thermal(assembly,n_VMs,varargin)
%  Author: Alexander Saccani, pHD candidate ETHZ, 03/2023

%this functions computes the VMs for the structure whose model must be
%provided in assembly, for different nodal temperatures configurations (the
%ones that are stored in vector Tsampl). 
% in varargin put in order the temperature samples and the temperature
% gradient (if required, e.g. in shell analysis)

T_samples = varargin{1};

if nargin > 3
 gradT_samples = varargin{2};
end

n_configs = size(T_samples,2);

nDofsF = assembly.Mesh.nDOFs; 
 
F = zeros(nDofsF,1);
M = mass_matrix(assembly);
MC = assembly.constrain_matrix(M);

%initialize output
VM = cell(n_configs,1);
static_sol = zeros(nDofsF,n_configs);
omega = zeros(n_VMs,n_configs);

for ii = 1:n_configs
    
    %temperature distribution T_ii
    T_ii = T_samples(:,ii);
    
    argTanStiff{1}= T_ii; %variable arguments for tangent stiffness matrix
    
    if nargin > 3
     gradT_ii = gradT_samples(:,ii);
     argTanStiff{2}= gradT_ii;
    end

    %compute static equilibrium
    [~,F0ii] = assembly.tangent_stiffness_and_force(zeros(nDofsF,1),argTanStiff{:}); %optimizable
    assembly.DATA.K = assembly.DATA.Kc;    % use cold stiffness to compute initial guess (better for convergence of nonlinear solution)
    assembly.DATA.F0 = F0ii;
    [~,eq_ii] = static_eq(assembly, F,'vararginTanStiffForce', {T_ii, gradT_ii}); %F0 and Kc are used inside this function to compute the initial guess

    static_sol(:,ii) = eq_ii; %static equilibrium

    %compute tangent stiffness matrix at static equilibrium (linearized model)
    [K_ii,~] = assembly.tangent_stiffness_and_force(eq_ii,argTanStiff{:});
    KC_ii = assembly.constrain_matrix(K_ii);
   

    %compute Vibration Modes for T_ii 
    [VM_ii,omega2_ii] = eigs(KC_ii,MC,n_VMs,'SM'); 
    omega(:,ii) = sqrt(diag(omega2_ii));

%     for jj = 1:n_VMs
%         VM_ii(:,jj) = VM_ii(:,jj)/max(sqrt(sum(VM_ii(:,jj).^2,2))); %normalize with respect to max. amplitude
%     end

    VM_ii = VM_ii./vecnorm(VM_ii); %normalize eigenvectors (they are unitary vectors)
    
    VM_ii = assembly.unconstrain_vector(VM_ii); %vibration modes for T_ii

    VM{ii} = VM_ii;
end

%reorder basis to avoid mode veering
% V_ref = VM{1,1};

% for ii = 2:n_configs
%     
%     P_ii = (VM{ii})'*V_ref;
%     [L,~,R] = svd(P_ii);
%     Q_ii = L*R';
%     VM{ii} = VM{ii}*Q_ii;
%     
% end

end