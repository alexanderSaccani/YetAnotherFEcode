classdef QUAD8ShellNeuHT < HTShellElementNeumann
    %BRICK 8 nodes brick element for heat transfer problem
    
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nNodes = 3          % number of nodes per element
        nDim = 1            % number of dimensions in local coordinates
        nelDOFs = 1
        nDOFPerNode = 1
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType = 'QUADNEU'
        thickness
    end
    
    
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = QUAD8ShellNeuHT(Material, thickness, Ngauss)
   
            if nargin == 2
                Ngauss = 2;
            end
            self.thickness = thickness;
            ShellElementConstructor(self, Material, Ngauss);
        end
        
        
        function J = jacobian_matrix(self, X)
            %______________________________________________________________
            %
            % Input: X is the vector of natural coordinates
            % Output: J is the jacobian matrix of physical coordinates (dim
            % = 3) with respect to the natural coordinates X (dim = 1).
            %          J = [dxdX;
            %               dydX;
            %               dzdX]
            %______________________________________________________________
           
            xyz = self.nodes; %3x3 matrix with xyz coordinates of nodes (stored row-wise)
            % shape function derivatives in natural coordinates
            dHn = 1/2*[ -1+2*X;   
                        1+2*X;
                        -4*X]';
            J = dHn*xyz;
            J = J.'; %Jacobian of transformation J = [dxdX
                     %                                dydX
                     %                                dzdX]


        end
        
        
        
    end
        
    
   
        
    methods (Static)
        
        function N = shape_functions(X)
            

            N = 1/2*[-X+X^2; X+X^2; 2-2*X^2];
              

        end
        
        function X = natural_coordinates
            X = [-1     % node 1
                  1      % node 2
                  0];	  % node 3

        end
        
    end
       
        
    
    
    
end

