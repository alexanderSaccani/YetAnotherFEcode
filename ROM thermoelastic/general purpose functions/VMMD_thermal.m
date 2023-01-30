function out = VMMD_thermal(assembly,T_samples,eq,n_VMs,logic_MD)

%this functions computes the VMs and MDs (if required) for the structure whose model must be
%provided in assebmby, for different nodal temperatures configurations (the
%ones that are stored in vector Tsampl). VMs can be computed around the
%undeformed configuration (if eq == 0 ) or with respect to the deformed by
%temperature static configuration (if eq == 1). In this case the output provides
%also the static deformation.

%set logic_MD to 0 if you don't want the modal derivatives in output.
%If you want them set logic_MD to 1

n_configs = size(T_samples,2);

nDofsF = assembly.Mesh.nDOFs; 
 
F = zeros(nDofsF,1);
M = mass_matrix(assembly);
MC = assembly.constrain_matrix(M);

%initialize output
VM = cell(n_configs,1);
static_eq = zeros(nDofsF,n_configs);
omega = zeros(n_VMs,n_configs);

if logic_MD == 1
    MD = cell(n_configs,1);
end

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
   
    %max(max(KC_ii-KC_ii.'))./max(max(KC_ii)) is tangent stiffness matrix
    %symmetric (check)
    
    %max(max(MC-MC.'))./max(max(MC));
    
    %compute Vibration Modes for T_ii 
    [VM_ii,omega2_ii] = eigs(KC_ii,MC,n_VMs,'SM'); 
    omega(:,ii) = sqrt(diag(omega2_ii));
   
    VM_ii = VM_ii./vecnorm(VM_ii); %normalize eigenvectors (they are unitary vectors)
    
    VM_ii = assembly.unconstrain_vector(VM_ii); %vibration modes for T_ii
    
    VM{ii} = VM_ii; %save output
    
    % compute modal derivatives if required for T_ii
    if logic_MD == 1
        
        %compute Modal Derivatives ii
        [MD_ii, names_ii] = modal_derivatives_thermal(assembly, assembly.Mesh.Elements, eq_ii, VM_ii, T_ii, 0);
        
        %form reduction basis
        MD_ii = assembly.constrain_vector(MD_ii);
        
        for jj = 1:size(MD_ii,2)
            MD_ii(:,jj) = MD_ii(:,jj)./vecnorm(MD_ii(:,jj)); %normalize eigenvectors (they are unitary vectors)
        end
        
        MD_ii = assembly.unconstrain_vector(MD_ii);
        
        MD{ii} = MD_ii; %save output
        
    end
    
    
end

%save output in struct format
out.VMs = VM;
out.omega = omega;
out.eq = static_eq;

if logic_MD == 1
    out.MDs = MD;
end

%reorder basis to avoid mode veering (this can be used to compare the space
%spanned by different basis corresponding to different temperature
%distributions
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