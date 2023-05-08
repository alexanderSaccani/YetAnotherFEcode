classdef Quad8Shell < ShallowShellElement
    %QUADSHELL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        
        nDOFPerNode = 5     % number of DOFs per node
        nNodes = 8          % number of nodes per element
        nDim = 2            % number of dimensions in local coordinates        
        nelDOFs
        
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        
        thickness
        
        membraneBending = []    %membrane bending relations (in the hypothesis of linear elastic material)      
        
        elType = 'QUAD'
        
        z = zeros(8,1); %define later, after the element table is created, assign height to shell elements with nodes from nodeIDs...
    end                 %if z is not set, it's default value is 0.
    
    
    methods
        
   % CONSTRUCTOR___________________________________________________________
        
        function obj = Quad8Shell(thickness, Material, varargin) %constructor
            
            if nargin < 3
                Ngauss = 2; %reduced integration
            else
                Ngauss = varargin{1};
            end
            
            obj.thickness = thickness;
            obj.nelDOFs = obj.nNodes * obj.nDOFPerNode;
            ShellElementConstructor(obj, Material, Ngauss);
            
            %initialize sparse matrices for extraction of membrane
            %displacements, rotations and vertical displacements from full vector of nodal displacements
            obj.initialize_extraction_matrices();
            
        end
        
  % SHAPE FUNCTIONS________________________________________________________
                
        function [fun,detJ,dNdx] = shape_function_derivatives(obj,X)
        
            xi = X(1); %natural coordinates
            eta = X(2);
            
            fun = 1/4*[...
                    -(1-xi)*(1-eta)*(1+xi+eta); 
                    -(1+xi)*(1-eta)*(1-xi+eta);
                    -(1+xi)*(1+eta)*(1-xi-eta); 
                    -(1-xi)*(1+eta)*(1+xi-eta);
                    2*(1-xi)*(1+xi)*(1-eta);  
                    2*(1-eta)*(1+eta)*(1+xi);
                    2*(1-xi)*(1+xi)*(1+eta);  
                    2*(1-eta)*(1+eta)*(1-xi)].'; 
            
            xy = obj.nodes; %xy coordinates of element's nodes
            
            dH = [ ...
                    -((eta - 1)*(xi + eta + 1))/4 - ((xi - 1)*(eta - 1))/4, -((xi - 1)*(xi + eta + 1))/4 - ((xi - 1)*(eta - 1))/4;
                    ((eta - 1)*(eta - xi + 1))/4 - ((xi + 1)*(eta - 1))/4,   ((xi + 1)*(eta - xi + 1))/4 + ((xi + 1)*(eta - 1))/4;
                    ((eta + 1)*(xi + eta - 1))/4 + ((xi + 1)*(eta + 1))/4,   ((xi + 1)*(xi + eta - 1))/4 + ((xi + 1)*(eta + 1))/4;
                    ((eta + 1)*(xi - eta + 1))/4 + ((xi - 1)*(eta + 1))/4,   ((xi - 1)*(xi - eta + 1))/4 - ((xi - 1)*(eta + 1))/4;
                                                        xi*(eta - 1),                                     xi^2/2 - 1/2;
                                                      1/2 - eta^2/2,                                      -eta*(xi + 1);
                                                       -xi*(eta + 1),                                     1/2 - xi^2/2;
                                                      eta^2/2 - 1/2,                                       eta*(xi - 1)]';

            J = dH*xy; %dx/dX 
            
            detJ = det(J); %determinant of the jacobian
            
            dNdx = J\dH;	% derivatives of shape function with respect to physical coordinates,
                      	% 2x16 matrix, [dNi_dx; dNi_dy]
                       	%  ( dg_dx = (dNdx*gnodes)' )               
        end
        
        %shape function for geometry description of midsurface 
        function fun = hxy(obj,X) 
            
            N = obj.shape_function(X);
                      
            fun = [N(1), 0, N(2), 0, N(3), 0, N(4), 0, N(5), 0, N(6), 0, N(7), 0, N(8), 0;...
                   0, N(1), 0, N(2), 0, N(3), 0, N(4), 0, N(5), 0, N(6), 0, N(7), 0, N(8)]; 
        end
        
        %height of shell midsurface function of natural coordinates X
        function out = hz(obj,X)     
            
            out = obj.shape_function(X)*obj.z; 
            
        end          

        %gauss to nodes matrices
        function out = gauss_to_nodes(obj)
            
            if (obj.quadrature.Ng == 2)
                
  out = [       0.1340   -0.5000   -0.5000    1.8660;
               -0.5000    1.8660    0.1340   -0.5000;
                1.8660   -0.5000   -0.5000    0.1340;
               -0.5000    0.1340    1.8660   -0.5000;
               -0.1830    0.6830   -0.1830    0.6830;
                0.6830    0.6830   -0.1830   -0.1830;
                0.6830   -0.1830    0.6830   -0.1830;
               -0.1830   -0.1830    0.6830    0.6830];

            else
                error('method to extract field at nodes from field at gauss points not specified for Ng > 2')
            end

        end
        
        %initialize extraction matrices
        function initialize_extraction_matrices(obj)
            
            EMu = sparse(8,40); %extraction matrix 
            EMu(1,1) = 1;       %this matrix, extracts the nodal u displacements from the vector of nodal displacements
            EMu(2,6) = 1;
            EMu(3,11) = 1;
            EMu(4,16) = 1;
            EMu(5,21) = 1;
            EMu(6,26) = 1;
            EMu(7,31) = 1;
            EMu(8,36) = 1;
            
            EMv = sparse(8,40);
            EMv(1,2) = 1;
            EMv(2,7) = 1;
            EMv(3,12) = 1;
            EMv(4,17) = 1;
            EMv(5,22) = 1;
            EMv(6,27) = 1;
            EMv(7,32) = 1;
            EMv(8,37) = 1;
            
            EMw = sparse(8,40);
            EMw(1,3) = 1;
            EMw(2,8) = 1;
            EMw(3,13) = 1;
            EMw(4,18) = 1;
            EMw(5,23) = 1;
            EMw(6,28) = 1;
            EMw(7,33) = 1;
            EMw(8,38) = 1;
            
            EMthx = sparse(8,40);
            EMthx(1,4) = 1;
            EMthx(2,9) = 1;
            EMthx(3,14) = 1;
            EMthx(4,19) = 1;
            EMthx(5,24) = 1;
            EMthx(6,29) = 1;
            EMthx(7,34) = 1;
            EMthx(8,39) = 1;
            
            EMthy = sparse(8,40);
            EMthy(1,5) = 1;
            EMthy(2,10) = 1;
            EMthy(3,15) = 1;
            EMthy(4,20) = 1;
            EMthy(5,25) = 1;
            EMthy(6,30) = 1;
            EMthy(7,35) = 1;
            EMthy(8,40) = 1;
            
            obj.initialization.EMu = EMu;
            obj.initialization.EMv = EMv;
            obj.initialization.EMw = EMw;
            obj.initialization.EMthx = EMthx;    
            obj.initialization.EMthy = EMthy;

            
        end

        
    end
    

       
    methods (Static)
        
        
        function fun = shape_function(X)
        
            xi = X(1); %natural coordinates
            eta = X(2);
            
            fun = 1/4*[...
                    -(1-xi)*(1-eta)*(1+xi+eta);  %shape fun 1D
                    -(1+xi)*(1-eta)*(1-xi+eta);
                    -(1+xi)*(1+eta)*(1-xi-eta); 
                    -(1-xi)*(1+eta)*(1+xi-eta);
                    2*(1-xi)*(1+xi)*(1-eta);  
                    2*(1-eta)*(1+eta)*(1+xi);
                    2*(1-xi)*(1+xi)*(1+eta);  
                    2*(1-eta)*(1+eta)*(1-xi)].';        
                  
        end
                     
                  
        end
        
         

        
             
        
end
    



