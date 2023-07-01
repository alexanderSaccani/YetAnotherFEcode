classdef HTContinuumElement < HeatTransferElement
    
    %continuum element for HEAT TRANSFER
    
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
    
    properties (Dependent) %WHAT ARE THESE FOR??
        area                % area of a 2D element (NaN otherwise)
        vol                	% volume of a 3D element (NaN otherwise)
    end
   
    methods (Abstract) % function that every continuum element must implement
        [G,detJ,dH] = shape_function_derivatives(self, X)
        N = shape_functions(X)
        X = natural_coordinates
    end
    
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
            
        end
        
        
        %Quadrature Rule
        function quadrature_rule(self, Ngauss)
            DIM = self.nDim;
            ELTYPE = self.elType;
            switch ELTYPE
                case {'HEX','QUAD'}
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
                                W(cont) = w(ii)*w(jj);
                                cont = cont+1;
                            end
                        end
                    end
                    self.quadrature.Ng = Ngauss;
                    self.quadrature.X = X;	% gauss integration points
                    self.quadrature.W = W;	% gauss integration weights
                case {'TET'}
                    [x, w] = inttet(Ngauss);
                    self.quadrature.Ng = Ngauss;
                    self.quadrature.X = x;	% gauss integration points
                    self.quadrature.W = w;	% gauss integration weights
                case {'WED'}
                    [w, x] = wedge_rule (Ngauss.lin, Ngauss.tri);
                    self.quadrature.Ng = Ngauss;
                    self.quadrature.X = x;	% gauss integration points
                    self.quadrature.W = w;	% gauss integration weights
                otherwise
                    error([' No quadrature rule found for the element "' ...
                        ELTYPE '". Accepted values are: HEX/QUAD, TET/TRI, WED'])
            end
        end
        
        
        
        function T = extract_element_data(self,T)
            
            %this function extracts temperatures at the element's nodes
            %from the vector of the nodal temperatures of the assembly
            T = T(self.nodeIDs);
            
        end
        
        
        %CONDUCTIVITY MATRIX_______________________________________________
        function K = conductivity_matrix(self)
            
            %get material properties
            k = self.Material.THERMAL_CONDUCTIVITY;
            
            %c = self.Material.SPECIFIC_HEAT;
            %rho = self.Material.DENSITY;
            
            %extract quadrature points and weights
            X = self.quadrature.X;
            W = self.quadrature.W;
            
            K = self.initialization.matrix;
            
            %gauss integration
            for ii = 1:length(self.quadrature.W)
                
                Xii = X(:,ii);
                Wii = W(ii);
                
                [detJ,dH] = shape_function_derivatives(self, Xii);
                
                dHx = dH(1,:); %[dN1dx,....,dNndx] 
                dHy = dH(2,:); %[dN1dy,....,dNndy]
                dHz = dH(3,:); %[dN1dz,....,dNndz]
                
                K = K + k*(dHx'*dHx + dHy'*dHy + dHz'*dHz)*detJ*Wii;
                    
            end

        end
        
        %THERMAL CAPACITY MATRIX___________________________________________
        function M = thermal_capacity_matrix(self)
            
            %extract material properties
            c = self.Material.SPECIFIC_HEAT;
            rho = self.Material.DENSITY;
            
            %extract quadrature points and weights
            X = self.quadrature.X;
            W = self.quadrature.W;
            
            M = self.initialization.matrix;
            
            %gauss integration
            for ii = 1:length(self.quadrature.W)
                
                Xii = X(:,ii);
                Wii = W(ii);
                
                [detJ,~] = shape_function_derivatives(self, Xii);
                
                N = self.shape_functions(Xii); % N=[N1,...Nn]
                
                M = M + c*rho*(N'*N)*detJ*Wii;
                    
            end
            
        end
        
        
        %LOAD MATRIX_______________________________________________________
        function LM = load_matrix(self)
            
            %this matrix multiplied with E allows to define the thermal load
            %coming from the internal generated heat per unit volume, E. E
            %is the vector of internal generated heat at element's nodes.
            
            %extract quadrature points and weights
            X = self.quadrature.X;
            W = self.quadrature.W;
            
            LM = self.initialization.matrix;
            
            %gauss integration
            for ii = 1:length(self.quadrature.W)
                
                Xii = X(:,ii);
                Wii = W(ii);
                
                [detJ,~] = shape_function_derivatives(self, Xii);
                
                N = self.shape_functions(Xii); % N=[N1,...Nn]
                
                LM = LM + (N'*N)*detJ*Wii;
                    
            end
                 
        end
        
        
    
    end
    
      
    
    
end

