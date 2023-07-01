classdef HTNeumannContinuumElement < HeatTransferElement
    
 %neumann continuum element for HEAT TRANSFER boundary conditions
    
    % ContinuumElement collects all the functions shared by all the
    % elements based on a continuum deformation model. It includes both
    % planar and solid elements, such as quadrilaterals, hexahedra, 
    % tetrahedra and wedges.
    
    properties (Abstract)
        nodes               % global coordinates of element nodes
        nodeIDs             % the index location of element nodes
        nNodes              % number of nodes per element
        nDim                % number of dimensions in local coordinates
        nelDOFs             % number of DOFs per element
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType              % 'HEX'/'QUAD'/'TET'/'TRI'/'WED'
                            % --> see method "quadrature_rule"!
    end
    
%     properties (Dependent) %WHAT ARE THESE FOR??
%         area                % area of a 2D element (NaN otherwise)
%         vol                	% volume of a 3D element (NaN otherwise)
%     end
   
%     methods (Abstract) % function that every continuum element must implement
%         [G,detJ,dH] = shape_function_derivatives(self, X)
%         N = shape_functions(X)
%         X = natural_coordinates
%     end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function ContinuumElementConstructor(self, Material, Ngauss)
            % _____________________________________________________________
            %
            % ContinuumElementConstructor(self, Material, Ngauss)
            % defines element's properties
            %______________________________________________________________
            self.Material = Material;
            ne = self.nelDOFs;
            nd = self.nDim;
            
            quadrature_rule(self, Ngauss);
            
            % INIZIALIZATION of some matrices (this should speedup
            % numerical integration)
            self.initialization.matrix = zeros(self.nelDOFs);
            self.initialization.vector = zeros(self.nelDOFs,1);
            
        end
        
        
        %Quadrature Rule
        function quadrature_rule(self, Ngauss)
            DIM = self.nDim;
            ELTYPE = self.elType;
            switch ELTYPE
                case {'HEXNEU','QUADNEU'} %Neumann Element
                    [x,w]=lgwt(Ngauss,-1,1);
                    X = zeros(DIM, Ngauss^DIM);
                    W = zeros(Ngauss^DIM, 1);
                    cont = 1;
                    for ii = 1:Ngauss
                        for jj = 1:Ngauss
                            if DIM == 3
                                for kk = 1:Ngauss
                                    X(:,cont) = [x(ii) x(jj) x(kk)].';
                                    W(cont) = w(ii)*w(jj)*w(kk);
                                    cont = cont+1;
                                end
                            elseif DIM == 2
                                X(:,cont) = [x(ii) x(jj)].';
                                W(:,cont) = w(ii)*w(jj);
                                cont = cont+1;
                            end
                        end
                    end
                    self.quadrature.Ng = Ngauss;
                    self.quadrature.X = X;	% gauss integration points
                    self.quadrature.W = W;	% gauss integration weights
%                 case {'TET'}
%                     [x, w] = inttet(Ngauss);
%                     self.quadrature.Ng = Ngauss;
%                     self.quadrature.X = x;	% gauss integration points
%                     self.quadrature.W = w;	% gauss integration weights
%                 case {'WED'}
%                     [w, x] = wedge_rule (Ngauss.lin, Ngauss.tri);
%                     self.quadrature.Ng = Ngauss;
%                     self.quadrature.X = x;	% gauss integration points
%                     self.quadrature.W = w;	% gauss integration weights
                otherwise
                    error([' No quadrature rule found for the element "' ...
                        ELTYPE '". Accepted values are: HEXNEU/QUADNEU']);% TET/TRI, WED'])
            end
        end
        
        
        
        
%         function T = extract_element_data(self,T)
%             
%             %this function extracts temperatures at the element's nodes
%             %from the vector of the nodal temperatures of the assembly
%             T = T(self.nodeIDs);
%             
%         end
        
        
        
        function Lth = thermal_flux_matrix(self)
            
            
            %extract quadrature points and weights
            X = self.quadrature.X;
            W = self.quadrature.W;
            
            Lth = self.initialization.matrix;          

            %gauss integration
            for ii = 1:length(self.quadrature.W)
                
                Xii = X(:,ii);
                Wii = W(ii);
                
                J = jacobian_matrix(self, Xii); %jacobian matrix. J = dxdX 
                %where x are physical coordinates of surface (or line) and
                %X are local natural coordinates of surface (or line).
                %dim(x) = dim(X) + 1
                
                %define the integration coefficient that is needed to
                %perform the boundary integral
                if self.nDim == 2 %for dimension 2 (integral on a surface)
                    intCoeff = norm([0,       -J(3,1),     J(2,1);
                                    J(3,1),   0,         -J(1,1);
                                    -J(2,1),  J(1,1),       0   ]*J(:,2));
%                     intCoeff1 = det(J(1:2,1:2));
%                     aa = intCoeff-intCoeff1;
                else
                    intCoeff = norm(J); %for dimension 1 (integral on a line)
                end

                N = self.shape_functions(Xii); %column vector
                
                Lth = Lth + N*N'*intCoeff*Wii;
                    
            end
               
        end
       
        
    
    
    
    end
    
end

