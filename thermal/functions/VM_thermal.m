function [VM,static_eq,omega] = VM_thermal(assembly,T_samples,eq,n_VMs)

%this functions computes the VMs for the structure whose model must be
%provided in assebmby, for different nodal temperatures configurations (the
%ones that are stored in vector Tsampl). VMs can be computed around the
%undeformed configuration (if eq == 0 ) or with respect to the deformed by
%temperature static configuration (if eq == 1). In this case the output provides
%also the static deformation.

n_configs = size(T_samples,2);

nDofsF = assembly.Mesh.nDOFs; 
 
F = zeros(nDofsF,1);
M = mass_matrix(assembly);
MC = assembly.constrain_matrix(M);

%initialize output
VM = cell(n_configs,1);
static_eq = zeros(nDofsF,n_configs);
omega = zeros(n_VMs,n_configs);

for ii = 1:n_configs
    %temperature distribution T_ii
    T_ii = T_samples(:,ii);

    if eq == 0
        eq_ii = zeros(nDofsF,1);
    else
        %compute static equilibrium
        [K_ii,~] = assembly.tangent_stiffness_and_force(zeros(nDofsF,1),T_ii); %optimizable
        assembly.DATA.K = K_ii;    % this value of K is used only for first guess solution (T not considered in the linear guess!)
        [~,eq_ii] = static_equilibrium_thermal(assembly, F, T_ii); %not super efficient since it computes also the linear solution
    end

    static_eq(:,ii) = eq_ii; %static equilibrium

    %compute tangent stiffness matrix at static equilibrium (linearized model)
    [K_ii,~] = assembly.tangent_stiffness_and_force(eq_ii,T_ii);
    KC_ii = assembly.constrain_matrix(K_ii);
    
    isKsym = max(max(K_ii-K_ii.'))/max(max(K_ii))*100

    %compute Vibration Modes for T_ii 
    [VM_ii,omega2_ii] = eigs(KC_ii,MC,n_VMs,'SM'); 
    omega(:,ii) = sqrt(diag(omega2_ii));

%     for jj = 1:n_VMs
%         VM_ii(:,jj) = VM_ii(:,jj)/max(sqrt(sum(VM_ii(:,jj).^2,2))); %normalize with respect to max. amplitude
%     end

% VM_ii = VM_ii./vecnorm(VM_ii); %normalize eigenvectors (they are unitary vectors)
    
    VM_ii = assembly.unconstrain_vector(VM_ii); %vibration modes for T_ii

    VM{ii} = VM_ii;
end


end