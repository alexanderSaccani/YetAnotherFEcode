classdef HEX8ElementHT < HTContinuumElement
    %BRICK 8 nodes brick element for heat transfer problem
    
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nNodes = 8          % number of nodes per element
        nDim = 3            % number of dimensions in local coordinates
        nelDOFs = 1
        nDOFPerNode = 1
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType = 'HEX'
    end
    
    properties
        thickness = 1       % element thickness  
    end
    
    
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = HEX8ElementHT(Material, Ngauss)
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
        
        function [detJ,dH] = shape_function_derivatives(self, X)
            %______________________________________________________________
            % 
            % [detJ,dH] = shape_function_derivatives(self,g,h,r)
            % INPUT:
            % self: element, X: natural coordinates
            % OUTPUT:
            % detJ = det(dxdxi) with x physical coordinates and xi natural
            % coordinates 
            % dH = shape function derivatives in physical coordinates
            % dH = [dN1dx,dN2dx,....,dN7dx,dN8dx;
            %       dN1dy,dN2dy,....,dN7dy,dN8dy;
            %       dN1dz,dN2dz,....,dN7dz,dN8dz]
            % dH is a 3x8 matrix
            %______________________________________________________________
            g = X(1);
            h = X(2);
            r = X(3);
            xyz = self.nodes;
            % shape functions derivatives (ABAQUS ORDER)
            dHn = [ -((h - 1)*(r - 1))/8    -((g - 1)*(r - 1))/8    -((g - 1)*(h - 1))/8
                     ((h - 1)*(r - 1))/8     ((g + 1)*(r - 1))/8     ((g + 1)*(h - 1))/8
                    -((h + 1)*(r - 1))/8    -((g + 1)*(r - 1))/8    -((g + 1)*(h + 1))/8
                     ((h + 1)*(r - 1))/8     ((g - 1)*(r - 1))/8     ((g - 1)*(h + 1))/8
                     ((h - 1)*(r + 1))/8     ((g - 1)*(r + 1))/8     ((g - 1)*(h - 1))/8
                    -((h - 1)*(r + 1))/8    -((g + 1)*(r + 1))/8    -((g + 1)*(h - 1))/8
                     ((h + 1)*(r + 1))/8     ((g + 1)*(r + 1))/8     ((g + 1)*(h + 1))/8
                    -((h + 1)*(r + 1))/8    -((g - 1)*(r + 1))/8    -((g - 1)*(h + 1))/8].';
            J = dHn*xyz;
            J1 = [0 0 0; 0 0 0; 0 0 0];
            J1(1,1) = J(2,2)*J(3,3) - J(2,3)*J(3,2);
            J1(2,1) = J(2,3)*J(3,1) - J(2,1)*J(3,3);
            J1(3,1) = J(2,1)*J(3,2) - J(2,2)*J(3,1);
            J1(1,2) = J(1,3)*J(3,2) - J(1,2)*J(3,3);
            J1(2,2) = J(1,1)*J(3,3) - J(1,3)*J(3,1);
            J1(3,2) = J(1,2)*J(3,1) - J(1,1)*J(3,2);
            J1(1,3) = J(1,2)*J(2,3) - J(1,3)*J(2,2);
            J1(2,3) = J(1,3)*J(2,1) - J(1,1)*J(2,3);
            J1(3,3) = J(1,1)*J(2,2) - J(1,2)*J(2,1);
            detJ = J(1,1)*J1(1,1) + J(1,2)*J1(2,1) + J(1,3)*J1(3,1);
            J1 = J1/detJ;
            dH = J1*dHn;        % derivatives in physical coordinates,
            % 3x8 matrix, [dNi_dx; dNi_dy; dNi_dz]
            % with i = 1...8
 
        end
        
        
        
    end
        
    
    
        
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h,r), 1x8 matrix
            % SHAPE FUNCTIONS FOR A 8-NODED HEXAHEDRON
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.4 Solid isoparametric quadrilaterals and hexahedra
            g = X(1);
            h = X(2);
            r = X(3);            
            N = 1/8*[(1-g)*(1-h)*(1-r); (1+g)*(1-h)*(1-r);
                     (1+g)*(1+h)*(1-r); (1-g)*(1+h)*(1-r);
                     (1-g)*(1-h)*(1+r); (1+g)*(1-h)*(1+r);
                     (1+g)*(1+h)*(1+r); (1-g)*(1+h)*(1+r)].';
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

