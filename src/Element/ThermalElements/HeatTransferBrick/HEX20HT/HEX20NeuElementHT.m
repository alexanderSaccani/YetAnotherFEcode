classdef HEX20NeuElementHT < HTNeumannContinuumElement
    %BRICK 8 nodes brick element for heat transfer problem
    
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nNodes = 8          % number of nodes per element
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
        
        function self = HEX20NeuElementHT(Material, Ngauss)
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
            g = X(1);
            h = X(2);
            xyz = self.nodes;
            % shape function derivatives in natural coordinates
            dHn = [ ...
                    -((h - 1)*(g + h + 1))/4 - ((g - 1)*(h - 1))/4, -((g - 1)*(g + h + 1))/4 - ((g - 1)*(h - 1))/4;
                    ((h - 1)*(h - g + 1))/4 - ((g + 1)*(h - 1))/4,   ((g + 1)*(h - g + 1))/4 + ((g + 1)*(h - 1))/4;
                    ((h + 1)*(g + h - 1))/4 + ((g + 1)*(h + 1))/4,   ((g + 1)*(g + h - 1))/4 + ((g + 1)*(h + 1))/4;
                    ((h + 1)*(g - h + 1))/4 + ((g - 1)*(h + 1))/4,   ((g - 1)*(g - h + 1))/4 - ((g - 1)*(h + 1))/4;
                                                        g*(h - 1),                                     g^2/2 - 1/2;
                                                      1/2 - h^2/2,                                      -h*(g + 1);
                                                       -g*(h + 1),                                     1/2 - g^2/2;
                                                      h^2/2 - 1/2,                                       h*(g - 1)]';
            J = dHn*xyz;
            J = J.';

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
            N = 1/4*[...
                    -(1-g)*(1-h)*(1+g+h); 
                    -(1+g)*(1-h)*(1-g+h);
                    -(1+g)*(1+h)*(1-g-h); 
                    -(1-g)*(1+h)*(1+g-h);
                    2*(1-g)*(1+g)*(1-h);  
                    2*(1-h)*(1+h)*(1+g);
                    2*(1-g)*(1+g)*(1+h);  
                    2*(1-h)*(1+h)*(1-g)];
        end
        
        function X = natural_coordinates
            X = [ ...
                -1  -1  % node 1 (corner)
                1   -1  % node 2 (corner)
                1   1	% node 3 (corner)
                -1  1   % node 4 (corner)
                0   -1  % node 5
                1   0   % node 6
                0   1   % node 7
                -1  0]; % node 8
        end
        
    end
       
        
    
    
    
end

