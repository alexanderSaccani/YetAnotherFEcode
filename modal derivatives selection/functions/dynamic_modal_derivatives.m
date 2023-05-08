% modal_derivatives
%
% Synthax:
% [MD, names] = modal_derivatives(myAssembly, elements, Phi, USEJULIA)
%
% Description: compute Modal Derivatives.
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode. It MUST contain the field
%     ".DATA.K" storing the unconstrained linear stiffness matrix K0
%   - elements: table of the elements
%   - Phi: tall matrix containing by columns the selected Vibration Modes
%     (unconstrained)
%   - USEJULIA (optional): if set to 1, the Julia module "DpROM.jl" is used
%     (default value is 0).
% OUTPUTS
%   - MD: tall matrix containing by columns all the Modal Derivatives that
%     can be computed from the given Phi (unconstrained MDs are returned).
%     MD(:,:,1) are the one corresponding to sum of frquencies, MD(:,:,2)
%     are the one corresponding to difference of frequencies
%   - names: matrix containing the subscripts of the MDs, for convenience.
%
%
% Created: Apr 2023
% Author: Alexander Saccani, Eignenossiche Technische Hochschule Zurich


function [MDs, names, angExcFreq] = dynamic_modal_derivatives(myAssembly, elements, VMs, omega)


n = size(VMs,1);
n_VMs = size(VMs,2);

K0 = myAssembly.DATA.K;
K0 = myAssembly.constrain_matrix( K0 );
M = myAssembly.DATA.M;
M = myAssembly.constrain_matrix( M );

MDs = zeros(n, n_VMs*(n_VMs+1)/2,2); %extra dimension is for frequency (sum of freq or diff)

names = zeros(n_VMs*(n_VMs+1)/2, 2);
angExcFreq = zeros(size(names,1),2);

kk = 1;
for jj = 1 : n_VMs
    
    Phi_j = VMs(:, jj);  
    
    t0 = tic;
    fprintf(' dKdq, assembling %d elements ...', size(elements,1))
    dK_deta_j = myAssembly.matrix('stiffness_derivative',Phi_j);
    fprintf(' %.2f s\n',toc(t0))
    
    dK_deta_j = myAssembly.constrain_matrix( dK_deta_j );
    
    %sum of frequencies
    for ii = 1 : n_VMs
        if ii < jj
            continue
        end
        
        omegai = omega(ii);
        omegaj = omega(jj);
        
        Phi_i = myAssembly.constrain_vector( VMs(:, ii) );
        
        %with sum of frequencies
        dPhi_i_deta_j = -1/2*(K0 -(omegai+omegaj)^2*M)\(dK_deta_j * Phi_i); 
        
%         a = find(K0);
%         num = K0(a)-(omegai+omegaj)^2*M(a);
%         KK = num./K0(a)*100;
%         
%         
%         KM = K0-(omegai+omegaj)^2*M;
%         diffKMK0 = KM-K0;
%         maxdiff = max(diffKMK0);
        
        th =  dPhi_i_deta_j / norm(dPhi_i_deta_j);
        MDs(:,kk,1) = myAssembly.unconstrain_vector( th );
        
        %with difference of frequencies
        dPhi_i_deta_j = -1/2*(K0 -(omegai-omegaj)^2*M)\(dK_deta_j * Phi_i); 
        
        th =  dPhi_i_deta_j / norm(dPhi_i_deta_j);
        MDs(:,kk,2) = myAssembly.unconstrain_vector( th );
        
        %names
        names(kk, :) = [ii jj];
        
        %excFreq 
        angExcFreq(kk,:) = [omegai+omegaj,abs(omegai-omegaj)];
        
        kk = kk + 1;
    end
    
    
end
