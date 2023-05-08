classdef ThermQuad4Shell < ThermShallowShellElement
    %QUADSHELL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        
        nDOFPerNode = 5     % number of DOFs per node
        nNodes = 4          % number of nodes per element
        nDim = 2            % number of dimensions in local coordinates        
        nelDOFs
        
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        
        thickness
        
        membraneBending = []    %membrane bending relations (in the hypothesis of linear elastic material)      
        
        elType = 'QUAD'
        
        z = zeros(4,1); %define later, after the element table is created, assign height to shell elements with nodes from nodeIDs...
    end                 %if z is not set, it's default value is 0.
    
    
    methods
        
   % CONSTRUCTOR___________________________________________________________
        
        function obj = ThermQuad4Shell(thickness, Material, varargin) %constructor
            
            if nargin < 3
                Ngauss = 2;
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
            
            fun = 1/4*[ +(xi - 1)*(eta - 1); %shape function (1d)
                      -(xi + 1)*(eta - 1);
                      +(xi + 1)*(eta + 1); 
                      -(xi - 1)*(eta + 1)].';    
            
            xy = obj.nodes; %xy coordinates of element's nodes
            
            dH = 1/4*[ eta-1,  xi-1; 
                        1-eta, -xi-1; 
                        eta+1,  xi+1; 
                       -eta-1,  1-xi].'; 
                   
            J = dH*xy; %dx/dX 
            
            detJ = det(J); %determinant of the jacobian
            
            dNdx = J\dH;	% derivatives of shape function with respect to physical coordinates,
                      	% 2x8 matrix, [dNi_dx; dNi_dy]
                       	% with i = 1...8 ( dg_dx = (dNdx*gnodes)' )               
        end
        
        %shape function for geometry description of midsurface 
        function fun = hxy(obj,X) 
            
            N = obj.shape_function(X);
            
            N1 = N(1);
            N2 = N(2);
            N3 = N(3);
            N4 = N(4);
            
            fun = [N1, 0, N2, 0, N3, 0, N4, 0;
                    0, N1, 0, N2, 0, N3, 0, N4]; 
        end
        
        %height of shell midsurface function of natural coordinates X
        function out = hz(obj,X)     
            
            out = obj.shape_function(X)*obj.z; 
            
        end
        
        
        %initialize extraction matrices
        function initialize_extraction_matrices(obj)
            
            EMu = sparse(4,20); %extraction matrix
            EMu(1,1) = 1;
            EMu(2,6) = 1;
            EMu(3,11) = 1;
            EMu(4,16) = 1;
            
            EMv = sparse(4,20);
            EMv(1,2) = 1;
            EMv(2,7) = 1;
            EMv(3,12) = 1;
            EMv(4,17) = 1;
            
            EMw = sparse(4,20);
            EMw(1,3) = 1;
            EMw(2,8) = 1;
            EMw(3,13) = 1;
            EMw(4,18) = 1;
            
            EMthx = sparse(4,20);
            EMthx(1,4) = 1;
            EMthx(2,9) = 1;
            EMthx(3,14) = 1;
            EMthx(4,19) = 1;
            
            EMthy = sparse(4,20);
            EMthy(1,5) = 1;
            EMthy(2,10) = 1;
            EMthy(3,15) = 1;
            EMthy(4,20) = 1;
            
            obj.initialization.EMu = EMu;
            obj.initialization.EMv = EMv;
            obj.initialization.EMw = EMw;
            obj.initialization.EMthx = EMthx;    
            obj.initialization.EMthy = EMthy;

            
        end

        %gauss to nodes matrices
        function out = gauss_to_nodes(obj)
            
            if (obj.quadrature.Ng == 2)
                
      out = [0.1340   -0.5000   -0.5000    1.8660;
            -0.5000    1.8660    0.1340   -0.5000;
            1.8660   -0.5000   -0.5000     0.1340;
            -0.5000    0.1340    1.8660   -0.5000];

            else
                error('method to extract field at nodes from field at gauss points not specified for Ng > 2')
            end

        end
             

        
    end    

       
    methods (Static)
        
        
        function fun = shape_function(X)
        
            xi = X(1); %natural coordinates
            eta = X(2);
            
            fun = 1/4*[ +(xi - 1)*(eta - 1); %shape function (1d)
                      -(xi + 1)*(eta - 1);
                      +(xi + 1)*(eta + 1); 
                      -(xi - 1)*(eta + 1)].';            
                  
        end        
                  
        end
            
        
end
    



