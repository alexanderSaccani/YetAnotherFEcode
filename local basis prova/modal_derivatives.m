% modal_derivatives

% Author: modified from Jacopo Marconi's code, Politecnico di Milano


function [MD, names] = modal_derivatives(myAssembly, Phi, K0)

n = size(Phi,1);
n_VMs = size(Phi,2);
eq = zeros(size(Phi,1),1);
eq = myAssembly.unconstrain_vector(eq);

MD = zeros(n, n_VMs*(n_VMs+1)/2);
names = zeros(n_VMs*(n_VMs+1)/2, 2);
kk = 1;
for jj = 1 : n_VMs
    
    Phi_j = myAssembly.unconstrain_vector(Phi(:, jj));

    dK_deta_j = myAssembly.matrix('stiffness_derivative',eq,Phi_j);
    
    dK_deta_j = myAssembly.constrain_matrix( dK_deta_j );
    
    for ii = 1 : n_VMs
        if ii < jj
            continue
        end
        
        Phi_i = Phi(:, ii);
        
        dPhi_i_deta_j = -K0\(dK_deta_j * Phi_i); 
        th =  dPhi_i_deta_j / max(abs(dPhi_i_deta_j));
        MD(:,kk) =  th ;
        names(kk, :) = [ii jj];
        kk = kk + 1;
    end
end
