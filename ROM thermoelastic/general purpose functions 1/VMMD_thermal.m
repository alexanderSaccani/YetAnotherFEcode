function out = VMMD_thermal(assembly,n_VMs,logic_MD,varargin)

%this functions computes the VMs and MDs (if required) for the structure whose model must be
%provided in assebmby, for different nodal temperatures configurations (the
%ones that are stored in vector Tsampl).

%set logic_MD to 0 if you don't want the modal derivatives in output.
%If you want them set logic_MD to 1

%vararginTanStiff is a cell vector with additional inputs (from displ) to
%compute the tangent stiffness matrix

%vararginStiffDer is a cell vector with additional inputs to compute the
%tangent stiffness derivative

%provide in varargin the T_samples and the gradT_samples if required (shell
%analysis)

if nargin < 4
    error('need to specify temperature samples')
end

T_samples = varargin{1};

flagGrad = 0;
if nargin > 4
   gradT_samples = varargin{2};
   flagGrad = 1;
end

n_configs = size(T_samples,2);

nDofsF = assembly.Mesh.nDOFs; 
 
Fnull = zeros(nDofsF,1);
M = mass_matrix(assembly);
MC = assembly.constrain_matrix(M);

%initialize output
VM = cell(n_configs,1);
equilibrium_point = zeros(nDofsF,n_configs);
omega = zeros(n_VMs,n_configs);

if logic_MD == 1
    MD = cell(n_configs,1);
end

for ii = 1:n_configs
    
    %temperature distribution T_ii (and gradient of T if required 
    T_ii = T_samples(:,ii);
    
    if flagGrad == 1
     gradT_ii = gradT_samples(:,ii);
     argTanStiff = {T_ii,gradT_ii};
    else
     argTanStiff = {T_ii};
    end
 

    %compute static equilibrium
    [~,F0ii] = assembly.tangent_stiffness_and_force(Fnull,argTanStiff{:}); %compute the constant stress due to T
    assembly.DATA.K = assembly.DATA.Kc;    % use cold stiffness to compute initial guess (better for convergence of nonlinear solution)
    assembly.DATA.F0 = F0ii;
    [~,eq_ii] = static_eq(assembly, Fnull, 'vararginTanStiffForce', argTanStiff); %F0 and Kc are used inside this function to compute the initial guess

    equilibrium_point(:,ii) = eq_ii; %static equilibrium

    %compute tangent stiffness matrix at static equilibrium (linearized model)
    [K_ii,~] = assembly.tangent_stiffness_and_force(eq_ii,argTanStiff{:});
    KC_ii = assembly.constrain_matrix(K_ii);
    
    %compute Vibration Modes for T_ii, gradT_ii
    [VM_ii,omega2_ii] = eigs(KC_ii,MC,n_VMs,'SM'); 
    omega(:,ii) = sqrt(diag(omega2_ii));
   
    VM_ii = VM_ii./vecnorm(VM_ii); %normalize eigenvectors (they are unitary vectors)
    
    VM_ii = assembly.unconstrain_vector(VM_ii); %vibration modes for T_ii
    
    VM{ii} = VM_ii; %save output
    
    % compute modal derivatives if required for T_ii, gradT_ii
    if logic_MD == 1
        
        %compute Modal Derivatives ii
        [MD_ii, names_ii] = modal_derivatives_thermal(assembly, assembly.Mesh.Elements, eq_ii, VM_ii, argTanStiff, {}, 0);
        
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
out.eq = equilibrium_point;

if logic_MD == 1
    out.MDs.data = MD;
    out.MDs.names = names_ii; 
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