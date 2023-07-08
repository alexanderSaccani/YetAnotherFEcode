classdef HEX8NeuElementHT < HTNeumannContinuumElement
    %BRICK 8 nodes brick element for heat transfer problem
    
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nNodes = 4          % number of nodes per element
        nDim = 2            % number of dimensions in local coordinates
        nelDOFs = 1
        nDOFPerNode = 1
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType = 'HEXNEU'
    end
    
    properties
        thickness = 1       % element thickness  
    end
    
    
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = HEX8NeuElementHT(Material, Ngauss)
            % _____________________________________________________________
            %
            % SELF-FUNCTION
            % self = Hex8Element(Material,Ngauss)
            % defines element's properties
            %______________________________________________________________
            if nargin == 1
                Ngauss = 2;
            end
            self.thickness = 1;
            ContinuumElementConstructor(self, Material, Ngauss);
        end
        
        
        function J = jacobian_matrix(self, X)
            %______________________________________________________________
            %
            % Input: X is the vector of natural coordinates
            % Output: J is the jacobian matrix of physical coordinates (dim
            % = 3) with respect to the natural coordinates (dim = 2).
            %          J = [dxdX1, dxdX2;
            %               dydX1, dydX2;
            %               dzdX1, dzdX2]
            %______________________________________________________________
            r = X(1);
            s = X(2);
            xyz = self.nodes; %4x3 matrix with xyz coordinates of nodes (stored row-wise)
            % shape function derivatives in natural coordinates
            dHn = 1/4*[ s-1,  r-1; 
                        1-s, -r-1; 
                        s+1,  r+1; 
                       -s-1,  1-r].';
            J = dHn*xyz;
            J = J.'; %Jacobian of transformation J = [dxdX1, dxdX2;
                     %                                dydX1, dydX2;
                     %                                dzdX1, dzdX2]
            
%             detJ = det(J);
%             dH = J\dHn;	% derivatives in physical coordinates,
%                       	% 2x8 matrix, [dNi_dx; dNi_dy]
%                        	% with i = 1...8
%             G = self.initialization.G;
%             G(1:2,1:2:end) = dH;
%             G(3:4,2:2:end) = dH;


        end
        
        
        
    end
        
    
   
        
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 8-NODED QUADRILATERAL
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.4 Solid isoparametric quadrilaterals and hexahedra
            g = X(1);
            h = X(2);
            N = 1/4*[ +(g - 1)*(h - 1); 
                      -(g + 1)*(h - 1);
                      +(g + 1)*(h + 1); 
                      -(g - 1)*(h + 1)];
        end
        
        function X = natural_coordinates
            X = [-1  -1   % node 1
                  1  -1   % node 2
                  1   1	  % node 3
                 -1   1]; % node 4
        end
        
    end
       
        
    
    
    
end

