classdef HTShellElement < HeatTransferElement
    
    % shell element for HEAT TRANSFER
    
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
                            % --> see method "quadrature_rule"
        thickness
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
        
        function ShellElementConstructor(self, Material, Ngauss)
            % _____________________________________________________________
            %
            % ContinuumElementConstructor(self, Material, Ngauss)
            % defines element's properties
            %______________________________________________________________
            self.Material = Material;
            
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
            
            %extract thickness
            t = self.thickness;
            
            %extract quadrature points and weights
            X = self.quadrature.X;
            W = self.quadrature.W;
            
            K = self.initialization.matrix;
            
            %gauss integration
            for ii = 1:length(self.quadrature.W)
                
                Xii = X(:,ii);
                Wii = W(ii);
                
                [dH,J] = shape_function_derivatives(self, Xii);
                J = J'; %3x2 matrix
                
                dHx = dH(1,:); %[dN1dx,....,dNndx] 
                dHy = dH(2,:); %[dN1dy,....,dNndy]
                dHz = dH(3,:); %[dN1dz,....,dNndz]
                
                %J is jacobian matrix. J = dxdX 
                %where x are physical coordinates of surface and
                %X are local natural coordinates of surface.
                %dim(x) = dim(X) + 1
                
                %define the integration coefficient that is needed to
                %perform the boundary integral

                intCoeff = norm([0,       -J(3,1),     J(2,1);
                                 J(3,1),   0,         -J(1,1);
                                 -J(2,1),  J(1,1),       0   ]*J(:,2));
                
                K = K + t*k*(dHx'*dHx + dHy'*dHy + dHz'*dHz)*intCoeff*Wii;
                    
            end

        end
        
        %THERMAL CAPACITY MATRIX___________________________________________
        function M = thermal_capacity_matrix(self)
            
            %extract material properties
            c = self.Material.SPECIFIC_HEAT;
            rho = self.Material.DENSITY;
            
            %extract thickness
            t = self.thickness;
            
            %extract quadrature points and weights
            X = self.quadrature.X;
            W = self.quadrature.W;
            
            M = self.initialization.matrix;
            
            %gauss integration
            for ii = 1:length(self.quadrature.W)
                
                Xii = X(:,ii);
                Wii = W(ii);
                
                [~,J] = shape_function_derivatives(self, Xii);
                J = J';
                
                N = self.shape_functions(Xii); % N=[N1,...Nn]
                
                intCoeff = norm([0,       -J(3,1),     J(2,1);
                    J(3,1),   0,         -J(1,1);
                    -J(2,1),  J(1,1),       0   ]*J(:,2));
                
                
                M = M + t*c*rho*(N'*N)*intCoeff*Wii;
                    
            end
            
        end
        
        
        %LOAD MATRIX_______________________________________________________
        function LM = load_matrix(self)
            
            %this matrix multiplied with Q allows to define the thermal load
            %coming from the heat flux [W/m^2] imposed to the structure per unit surface 
            % (the sum of the heat from upper and lower surface of the shell)
            %Q is the vector of heat transfered to the shell through its boundry at the element's nodes.
            
            %extract quadrature points and weights
            X = self.quadrature.X;
            W = self.quadrature.W;
            
            LM = self.initialization.matrix;
            
            
            %gauss integration
            for ii = 1:length(self.quadrature.W)
                
                Xii = X(:,ii);
                Wii = W(ii);
                
                [~,J] = shape_function_derivatives(self, Xii);
                J = J';
                
                N = self.shape_functions(Xii); % N=[N1,...Nn]
                
                intCoeff = norm([0,       -J(3,1),     J(2,1);
                    J(3,1),   0,         -J(1,1);
                    -J(2,1),  J(1,1),       0   ]*J(:,2));
                
                detJ = det(J(1:2,1:2));

                LM = LM + (N'*N)*intCoeff*Wii;
                    
            end
                 
        end
        
        
    
    end
    
      
    
    
end

