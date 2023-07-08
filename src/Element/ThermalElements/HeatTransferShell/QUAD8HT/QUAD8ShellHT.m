classdef QUAD8ShellHT < HTShellElement
    %Shell 8 nodes element for heat transfer problem
    
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
        elType = 'QUAD'
        
    end
    
    properties
        thickness = 1       % element thickness  
    end
    
    
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = QUAD8ShellHT(Material, thickness, Ngauss)
            % _____________________________________________________________
            %
            % SELF-FUNCTION
            % self = QUAD8ShellHT(Material,thickness,Ngauss)
            % defines element's properties
            %______________________________________________________________
            if nargin == 2
                Ngauss = 2;
            end
            self.thickness = thickness;
            ShellElementConstructor(self, Material, Ngauss);
        end
        
       
        
        function [dH,J] = shape_function_derivatives(self, X)
            %______________________________________________________________
            %
            % [dH,J] = shape_function_derivatives(self, X)
            %
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
            J = dHn*xyz; %2x3 matrix, dxdX'
            dH = J\dHn;	% derivatives in physical coordinates,
                      	% 3x8 matrix, [dNi_dx; dNi_dy; dNi_dz]
                       	% with i = 1...8

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
                    2*(1-h)*(1+h)*(1-g)]';
        end
        
        function X = natural_coordinates
            X = [-1 -1 -1
                 +1 -1 -1
                 +1 +1 -1
                 -1 +1 -1
                 -1 -1 +1
                 +1 -1 +1
                 +1 +1 +1
                 -1 +1 +1];
        end
        
    end
       
        
    
    
    
end

