function out = decodeDofsNodes(D,nNodes,nDOFPerNode)
%DECODEDOFSNODES this function partitions a snapshot matrix (2d) of displacements /
%forces in a 3d array in which it is easy to access with # node and # of
%dof for specific node to given element. ATTENTION: It works if the same elements is
%used for discretization.

% INPUT: D is snapshot matrix of displacements / forces. Snapshots are stored
% in columns
%        nNodes
%        nDOFperNode
%
% OUTPUT: out(ii,jj,kk) ii-># node, jj-> # dof (for given node, kk-> #
% snapshot

n_tsamp = size(D,2);

out = zeros(nNodes,nDOFPerNode,n_tsamp); % (node, dof of node, sample)

for ii = 1:n_tsamp
   disp_ii = D(:,ii);
   out(:,:,ii) = (reshape(disp_ii,nDOFPerNode,[])).';
end

end

