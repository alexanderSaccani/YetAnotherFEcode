classdef HTShellElementNeumann< HeatTransferElement
    
 %neumann continuum element for HEAT TRANSFER boundary conditions
    
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
        thickness
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
                        if (DIM == 3 || DIM == 2)
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
                        elseif DIM == 1
                            X(:,cont) = x(ii);
                            W(cont) = w(ii);
                            cont = cont+1;
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
            
            t = self.thickness;
            
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

                intCoeff = norm(J); %for dimension 1 (integral on a line)

                N = self.shape_functions(Xii); %column vector
                
                Lth = Lth + t*(N*N')*intCoeff*Wii;
                    
            end
               
        end
       
        
    
    
    
    end
    
end

