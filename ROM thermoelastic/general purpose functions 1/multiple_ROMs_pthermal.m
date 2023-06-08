function [ROMs,VMs,MDs,eqDisp,natFreq] = multiple_ROMs_pthermal(fullAssembly, p, Tfun, gradFun, n_VMs, thermEq, orthog, reordBasis, unconstrainBasis)
%  Author: Alexander Saccani, pHD candidate ETHZ, 12/2022
%  last modified: 25.05.2023
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
%   2) p is the parameter vector from which depends the temperature
%   configuration. Each column represent a new value of the parameter
%   vector
%   3) Tfun(p) is a fun(p(:,ii)) that returns with the samples of T
%   corresponding to a given parameter
%   4) gradFun. struct array with fields flag and fun. If the temperature
%   flag is true, the temperature gradient is also considered and is evaluated as function of gradFun.param. 
%   Pay attention: the temperature gradient is define only for shell.
%   5) n_VMs: number of VMs that are included in the ROBs along with the
%   corresponding MDs.
%   6) thermEq is a logical. If true the thermal eq is included in the
%   basis (first element)
%   7) orthog is a logical. If true the basis is orthogonalized
%   8) reordBasis is a logical. If true the basis is reorder to avoid mode
%   veering
%   9) unconstrainBasis is a logical. If true the basis is unconstrained
%
%  OUTPUT:
%   1)ROMs: cell array containg in each cell a ROM corresponding to
%   different temperature configs. Each cell has fields:
%               1) thermal_eq
%               2) omega
%               3) V
%   2)VMs: cell array containing in each cell the VMs corresponding to a
%   different thermal parameter (unconstrained vector)
%   3)MDs: cell array containing in each cell the MDs corresponding to a
%   different thermal parameter (unconstrained vector)

n_ROMs = size(p,2);

nDofsF = fullAssembly.Mesh.nDOFs; 


F = zeros(nDofsF,1);
M = mass_matrix(fullAssembly);
MC = fullAssembly.constrain_matrix(M);

%initialize output
ROMs = cell(n_ROMs,1);
VMs = cell(n_ROMs,1);
MDs = cell(n_ROMs,1);
eqDisp = cell(n_ROMs,1);
natFreq = zeros(n_VMs,n_ROMs);

%remove field F0 (you don't know for which temp. config. it was initialized!
if isfield(fullAssembly.DATA,'F0')
    fullAssembly.DATA = rmfield(fullAssembly.DATA,'F0'); 
end

for ii = 1:n_ROMs

    %temperature distribution T_ii
%     arginTfun = num2cell(p(:,ii));
%     T_ii= Tfun(arginTfun{:});
    T_ii= Tfun(p(:,ii));

    %compute static equilibrium
    [~,F_ii] = tang_stiff_force(zeros(nDofsF,1),T_ii,ii); %optimizable
    fullAssembly.DATA.F0 = F_ii;   %static force deriving from T distribution
    fullAssembly.DATA.K =  fullAssembly.DATA.Kc;   % this value of K is used only for first guess solution (Kc cold is better for convergence of static analysis) ;
    [~,u_eq_ii] = equilibrium_state(F, T_ii, ii); %not super efficient since it computes also the linear solution
    ROMs{ii}.thermal_eq = fullAssembly.constrain_vector(u_eq_ii); %save in output
    eqDisp{ii} = u_eq_ii; %save thermal eq. corresponding to configuration ii

    %compute tangent stiffness matrix at static equilibrium (linearized model)
    [K_ii,~] = tang_stiff_force(u_eq_ii,T_ii,ii);
    KC_ii = fullAssembly.constrain_matrix(K_ii);

    %compute Vibration Modes for T_ii 
    [VM_ii,omega2_ii] = eigs(KC_ii,MC,n_VMs,'SM'); 
    omega_ii = sqrt(diag(omega2_ii));
    for jj = 1:n_VMs
        VM_ii(:,jj) = VM_ii(:,jj)./vecnorm(VM_ii(:,jj)); %normalize eigenvectors (they are unitary vectors)
    end
    VM_ii = fullAssembly.unconstrain_vector(VM_ii); %vibration modes for T_ii
    VMs{ii} = VM_ii; %store VMs at temp config ii
    
    %compute Modal Derivatives ii
    [MD_ii, ~] = modal_der(u_eq_ii,T_ii,VM_ii,ii);
    for jj = 1:size(MD_ii,2)
        MD_ii(:,jj) = MD_ii(:,jj)./vecnorm(MD_ii(:,jj)); %normalize eigenvectors (they are unitary vectors)
    end
    MDs{ii} = MD_ii; %store MDs at temp config ii
    
    %form reduction basis
    MD_ii = fullAssembly.constrain_vector(MD_ii);
    VM_ii = fullAssembly.constrain_vector(VM_ii); 
    
    if thermEq
        V_ii = [fullAssembly.constrain_vector(u_eq_ii),VM_ii,MD_ii]; %include the thermal eq in the basis
    else
        V_ii = [VM_ii,MD_ii]; %form ROB
    end
    
    %uncostrain basis
    if unconstrainBasis
        V_ii = fullAssembly.unconstrain_vector(V_ii);
    end
    
    %store output (constrained ROB)
    if orthog
        ROMs{ii}.V = orth(V_ii); %orthgonalize basis if required
    else
        ROMs{ii}.V = V_ii;
    end
    
    
    ROMs{ii}.omega = omega_ii; %store values of angular frequencies
    natFreq(:,ii) = omega_ii/2/pi;

end

% reorder basis to avoid mode veering if required (make sure that orth is
% true)
if ( reordBasis && orthog )
    V_ref = ROMs{1}.V;
    for ii = 2:n_ROMs

        P_ii = (ROMs{ii}.V)'*V_ref;
        [L,~,R] = svd(P_ii);
        Q_ii = L*R';
        ROMs{ii}.V = (ROMs{ii}.V)*Q_ii;

    end
elseif reordBasis && ~orthog
    warning('In order to reoder basis set orthog flag to true')
end


%this function allows to pass the gradient to the tangent stiffness and
%force fucntion at element level (if required)
    function [Ktt,F] = tang_stiff_force(disp,T,nPsample)
        
        if gradFun.flag 
            
            gradTFun = gradFun.fun;
            pGrad = gradFun.param;
            gradT = gradTFun(pGrad(:,nPsample));
            [Ktt,F] = fullAssembly.tangent_stiffness_and_force(disp,T,gradT);
            
        else
            
            [Ktt,F] = fullAssembly.tangent_stiffness_and_force(disp,T);
            
        end
        
    end

%this function allows to pass the gradient for the static equilibrium analysis 
    function [ulin,unl] = equilibrium_state(Fext,T,nPsample)
        
        if gradFun.flag 
            
            gradTFun = gradFun.fun;
            pGrad = gradFun.param;
            gradT = gradTFun(pGrad(:,nPsample));
            %in future consider to put F0 as field of fullAssmebly to speed
            %up computations
            [ulin,unl] = static_eq( fullAssembly, Fext, 'vararginTanStiffForce', {T,gradT});
            
        else
            %in future consider to put F0 as field of fullAssmebly to speed
            %up computations
            [ulin,unl] = static_eq( fullAssembly, Fext, 'vararginTanStiffForce', {T});
            
        end
        
    end

%this function allows to compute the modal derivatives 
    function [MDs,namesMDs] = modal_der(dEq,T,VMs,nPsample)
        
        %dEq = equilibrium configuration of expansion
        if gradFun.flag 
            gradTFun = gradFun.fun;
            pGrad = gradFun.param;
            gradT = gradTFun(pGrad(:,nPsample));
            varargTanStiff = {T,gradT};
            varargTanStiffDer = {}; %pay attention for beam elements
            [MDs, namesMDs] = modal_derivatives_thermal(fullAssembly, fullAssembly.Mesh.Elements, dEq, VMs, varargTanStiff, varargTanStiffDer, 0);
        else
            varargTanStiff = {T};
            varargTanStiffDer = {};
            [MDs, namesMDs] = modal_derivatives_thermal(fullAssembly, fullAssembly.Mesh.Elements, dEq, VMs, varargTanStiff, varargTanStiffDer, 0);
            
        end
        
    end



end

