function ROMs = multiple_ROMs_thermal(fullAssembly, T_samples, n_VMs)
%  Author: Alexander Saccani, pHD candidate ETHZ, 12/2022
%
%  MULTIPLE_ROMS_THERMAL 
%  Output. Cell array containing ROMs for each T distribution in the
%  columns of T_samples. Each element of cell array has fields V and omega.
%  Starting from  the fullAssebly, nodal temperatures (stored in T_samples) 
%  corresponding to different temperature configurations, this function
%  computes multiple linearized ROMs with specified number of VMs and MDs
%  starting from the thermal equilibrium point. The basis are aligned
%  starting to a reference basis (basis for the first ROM) to allow 
%  interpolation.

%  INPUT:
%   1) fullAssembly
%   2) T_samples: array containing in each column nodal vector samples
%   corresponding to different temperature configurations
%   3) n_VMs: number of VMs that are included in the ROBs along with the
%   corresponding MDs.
%
%  OUTPUT:
%   1)ROMs: cell array containg in each cell a ROM corresponding to
%   different temperature configs. Each cell has fields:
%               1) thermal_eq
%               2) omega
%               3) V

n_ROMs = size(T_samples,2);

nDofsF = fullAssembly.Mesh.nDOFs; 


F = zeros(nDofsF,1);
M = mass_matrix(fullAssembly);
MC = fullAssembly.constrain_matrix(M);

%initialize output
ROMs = cell(n_ROMs,1);

for ii = 1:n_ROMs

    %temperature distribution T_ii
    T_ii= T_samples(:,ii);

    %compute static equilibrium
    [K_ii,~] = fullAssembly.tangent_stiffness_and_force(zeros(nDofsF,1),T_ii); %optimizable
    fullAssembly.DATA.K = K_ii;    % this value of K is used only for first guess solution (T not considered in the linear guess!)
    [~,u_eq_ii] = static_equilibrium_thermal(fullAssembly, F, T_ii); %not super efficient since it computes also the linear solution
    ROMs{ii}.thermal_eq = fullAssembly.constrain_vector(u_eq_ii); %save in output

    %compute tangent stiffness matrix at static equilibrium (linearized model)
    [K_ii,~] = fullAssembly.tangent_stiffness_and_force(u_eq_ii,T_ii);
    KC_ii = fullAssembly.constrain_matrix(K_ii);

    %compute Vibration Modes for T_ii 
    [VM_ii,omega2_ii] = eigs(KC_ii,MC,n_VMs,'SM'); 
    omega_ii = sqrt(diag(omega2_ii));
    for jj = 1:n_VMs
        %VM_ii(:,jj) = VM_ii(:,jj)/max(sqrt(sum(VM_ii(:,jj).^2,2))); %normalize with respect to max. amplitude
        VM_ii(:,jj) = VM_ii(:,jj)./vecnorm(VM_ii(:,jj)); %normalize eigenvectors (they are unitary vectors)
    end
    VM_ii = fullAssembly.unconstrain_vector(VM_ii); %vibration modes for T_ii

    %compute Modal Derivatives ii
    [MD_ii, ~] = modal_derivatives_thermal(fullAssembly, fullAssembly.Mesh.Elements, u_eq_ii, VM_ii, T_ii, 0);
    for jj = 1:size(MD_ii,2)
        MD_ii(:,jj) = MD_ii(:,jj)./vecnorm(MD_ii(:,jj)); %normalize eigenvectors (they are unitary vectors)
    end
    
    %form reduction basis
    MD_ii = fullAssembly.constrain_vector(MD_ii);
    VM_ii = fullAssembly.constrain_vector(VM_ii); 
    
    V_ii = [VM_ii,MD_ii]; %form ROB
    %V_ii = VM_ii; %uncomment if you want a basis with only MDs
    
    ROMs{ii}.V = orth(V_ii); %store output (constrained ROB)
    ROMs{ii}.omega = omega_ii; %store values of angular frequencies

    %reduced matrices (do in future to allow interpolation of matrices) instead
    %of basis)
    % ROMs(ii).Mr = V_ii'*M*V_ii;
    % ROMs(ii).Kr = V_ii'*K_ii*V_ii;
    % in future add also K2,K3 for each model

end

% reorder basis to avoid mode veering
V_ref = ROMs{1}.V;

for ii = 2:n_ROMs
    
    P_ii = (ROMs{ii}.V)'*V_ref;
    [L,~,R] = svd(P_ii);
    Q_ii = L*R';
    ROMs{ii}.V = (ROMs{ii}.V)*Q_ii;
    
end



end

