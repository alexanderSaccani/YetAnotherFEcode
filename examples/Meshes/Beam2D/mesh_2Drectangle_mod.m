% mesh_2Drectangle
%
% syntax:
% [nodes, elements, bset] = ...
%           mesh_2Drectangle(Lx,Ly,nx,ny)
%
% Description: create the mesh for a rectangle using quadratic
%              quadrilaterals (Quad8).
%
% INPUTS    Lx, Ly: length, width
%           nx,ny: number of elements for Lx and Ly, respectively.
%           
% OUTPUTS   nodes --> [x, y]
%           elements --> [node_1, ... , node_8]
%           bset --> struct containing node sets (IDs) of the 4 external
%           edges. They are ordered as: x=0, y=0, x=l, y=w.
%           belemens --> struct containing element connectivity of the
%           nodes in bset

function [nodes, elements, bset, belements] = mesh_2Drectangle_mod(Lx,Ly,nx,ny,elementType)

nel = nx*ny;
elementType = upper(elementType);

fprintf([' Meshing %d elements (' elementType ') ... '], nel)
tic

switch elementType
    case 'QUAD4'
        nNodes = 4;
        elnodes = [         % element coordinates in the natural space (/2)
         0     0
         1     0
         1     1
         0     1];
        
        %Element's boundary node IDs for boundary element construction
        nNodesBoundaryEl = 2; %number of nodes of each edge of the element
        IDsNodesBoundary{1} = [1,4];
        IDsNodesBoundary{2} = [1,2];
        IDsNodesBoundary{3} = [2,3];
        IDsNodesBoundary{4} = [3,4];
    case 'QUAD8'
        nNodes = 8;
        elnodes = [         % element coordinates in the natural space (/2)
             0     0
             1     0
             1     1
             0     1
             0.5   0
             1     0.5
             0.5   1
             0     0.5];
         %Element's boundary node IDs for boundary element construction
         nNodesBoundaryEl = 3;  %number of nodes of each edge of the element
         IDsNodesBoundary{1} = [1,4,8];
         IDsNodesBoundary{2} = [1,2,5];
         IDsNodesBoundary{3} = [2,3,6];
         IDsNodesBoundary{4} = [3,4,7];
end
lx = Lx/nx;
ly = Ly/ny;
elnodes(:,1) = elnodes(:,1)*lx; % element coordinates in physical space
elnodes(:,2) = elnodes(:,2)*ly;
% create nodes
nodes = zeros(nel*nNodes,2);
nn = 1;

%initialize quantities for boundary element set
nodeSet1 = zeros(nx,nNodesBoundaryEl);
nodeSet2 = zeros(ny,nNodesBoundaryEl);
nodeSet3 = nodeSet1;
nodeSet4 = nodeSet2;

countSet1 = 1;
countSet2 = 1;
countSet3 = 1;
countSet4 = 1;
        
for ii = 1:nx
    for jj = 1:ny   
        elnodes_temp = elnodes;
        elnodes_temp(:,1) = elnodes_temp(:,1)+lx*(ii-1);
        elnodes_temp(:,2) = elnodes_temp(:,2)+ly*(jj-1);
        nodes(nn:nn+nNodes-1,:)  = elnodes_temp;
        
        if ii == 1
            nodeSet1(countSet1,:) = nn-1+IDsNodesBoundary{1};
            countSet1 = countSet1+1;
        elseif ii == nx
            nodeSet3(countSet3,:) = nn-1+IDsNodesBoundary{3};
            countSet3 = countSet3+1;
        end
        
        if jj == 1
            nodeSet2(countSet2,:) = nn-1+IDsNodesBoundary{2};
            countSet2 = countSet2+1;
        elseif jj == ny
            nodeSet4(countSet4,:) = nn-1+IDsNodesBoundary{4};
            countSet4 = countSet4+1;
        end
        
        
        
        nn = nn+nNodes;
    end
end

tol = 1e12;
nodes = round(nodes*tol) / tol;

% remove duplicate nodes from 'nodes' matrix
[nodes, ~, ic] = unique(nodes, 'rows', 'stable');
elements = ic;
elements = reshape(elements, nNodes, nel)';
idb = reshape(1 : length(nodes(:)), 2, size(nodes,1))';

conn = zeros(nel, nNodes*2);
for ii = 1:nel
    conn(ii,:) = reshape(idb(elements(ii,:),:)',1, nNodes*2);
end

%extract the element connectivity of boundary elements
nodeSet1 = reshape(nodeSet1',[],1);
nodeSet1 = ic(nodeSet1);
elementSet1 = (reshape(nodeSet1,nNodesBoundaryEl,[]))';

nodeSet2 = reshape(nodeSet2',[],1);
nodeSet2 = ic(nodeSet2);
elementSet2 = (reshape(nodeSet2,nNodesBoundaryEl,[]))';

nodeSet3 = reshape(nodeSet3',[],1);
nodeSet3 = ic(nodeSet3);
elementSet3 = (reshape(nodeSet3,nNodesBoundaryEl,[]))';

nodeSet4 = reshape(nodeSet4',[],1);
nodeSet4 = ic(nodeSet4);
elementSet4 = (reshape(nodeSet4,nNodesBoundaryEl,[]))';


% external surfaces
tol = 1e-15;
node_IDs = 1 : size(nodes, 1);
edge1_X0_nodes = node_IDs( nodes(:,1)==0 );
edge2_Y0_nodes = node_IDs( nodes(:,2)==0 );
edge3_XL_nodes = node_IDs( abs(nodes(:,1)-Lx)<tol );
edge4_YW_nodes = node_IDs( abs(nodes(:,2)-Ly)<tol );

% edge1_X0_dofs=reshape(idb(edge1_X0_nodes,2:end)',1,length(edge1_X0_nodes)*2);
% edge2_Y0_dofs=reshape(idb(edge2_Y0_nodes,2:end)',1,length(edge2_Y0_nodes)*2);
% edge3_XL_dofs=reshape(idb(edge3_XL_nodes,2:end)',1,length(edge3_XL_nodes)*2);
% edge4_YW_dofs=reshape(idb(edge4_YW_nodes,2:end)',1,length(edge4_YW_nodes)*2);

bset = {edge1_X0_nodes, edge2_Y0_nodes, edge3_XL_nodes, edge4_YW_nodes};
belements = {elementSet1, elementSet2, elementSet3, elementSet4}; %element connectivity of boundary elements

fprintf(' %.2f s\n',toc)
fprintf('  Nodes: %d \n',size(nodes,1))
fprintf('  Dofs:  %d \n\n',size(nodes,1)*2)

