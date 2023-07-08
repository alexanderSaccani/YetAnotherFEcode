classdef ShallowShellElement < Element
    
    %SHELLELEMENT shell element superclass
    % Thermal shallow shell element with 5 dofs per node, following the
    % Mindlin-Reissner kinematic formulation for displacements and Von
    % Karman strain displacements relation
    
    properties (Abstract)
        nodes               % global coordinates of element nodes
        nodeIDs             % the index location of element nodes
        nDOFPerNode         % number of DOFs per node
        nNodes              % number of nodes per element
        nDim                % number of dimensions in local coordinates
        nelDOFs             % number of DOFs per element
        
        quadrature      	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType              % ''/'QUAD','TRI'
                            % --> see method "quadrature_rule"
        thickness           % thickness of element
        
        z                   % height of nodes   
    end
    
    methods
        
        function ShellElementConstructor(obj, Material, Ngauss)
         
                % ShellElementConstructor(self, Material, Ngauss)
                % defines element's properties
                
                obj.Material = Material;
       
                
                ne = obj.nelDOFs;
                
                %define the quadrature rule
                quadrature_rule(obj, Ngauss); 
                
                %computation of the properties through the thickness (in
                %future it is possible to use material that varies through
                %the section
                C = get_stress_strain_matrix_2D(obj.Material);
                
%                W = obj.quadrature.W; %weight for numerical 
%                X = obj.quadrature.X; %integration points

                %membrane bending constitutive law
                t = obj.thickness;
                obj.initialization.Cm = C*t;
                obj.initialization.Cb = C*t^3/12 ; 
                
                %intialize some useful quantities 
                obj.initialization.K = zeros(ne); 	 % stiffness-element matrix
                obj.initialization.M = zeros(ne);    % mass matrix
                obj.initialization.F = zeros(ne,1); % element internal force
                
                alpha = obj.Material.THERMAL_EXPANSION_COEFFICIENT;
                E = obj.Material.YOUNGS_MODULUS;
                nu = obj.Material.POISSONS_RATIO;
                Gshear = E/2/(1+nu); %linear elastic material
                shearFactor = 5/6; %from Crisfield (or pi^2/12)
                
                obj.initialization.alpha = [alpha;alpha;0]; %[alphax,alphay,0] the thermal exp. coeff does not contribute to shear
                obj.initialization.Cgamma = eye(2)*Gshear*obj.thickness*shearFactor;
                
                H = [1,0,0,0;0,0,0,1;0,1,1,0];
                obj.initialization.H = H;        	 % linear membrane strain = H*G_{uv}*p
                
                
                obj.initialization.AEfun = @(dwdx,dzdx) ...
                    [1/2*dwdx(1)+dzdx(1)   0                  ;
                     0                     1/2*dwdx(2)+dzdx(2);
                    dwdx(2)+dzdx(2)        dzdx(1)            ]; %nonlinear membrane strain = AEfun(dwdx,dzdx)*G_w*p
                
                obj.initialization.AdelEfun = @(dwdx,dzdx) ...
                    [dwdx(1)+dzdx(1)       0              ;
                     0                     dwdx(2)+dzdx(2);
                     dwdx(2)+dzdx(2)       dwdx(1)+dzdx(1)]; % nonlinear deltaE = AdelEfun(dwdx,dzdx)*Gw*dp
                 
                obj.initialization.zMatrixFun = @(dzdx)...
                    [dzdx(1)    0   ;
                      0      dzdx(2);
                     dzdx(2) dzdx(1)];
                 
                % dE = (Hprime + (lambdaDelE*p)*G.Gw)dp, 
                % Hprime = zMatrix*G.Gw + H*G.Guv
                LambdaDelE = tenzeros(3,2,2);
                LambdaDelE(1,1,1) = 1; LambdaDelE(2,2,2) = 1; 
                LambdaDelE(3,1,2) = 1; LambdaDelE(3,2,1) = 1;
                obj.initialization.LambdaDelE = LambdaDelE;
                
                % E = (Hprime + (lambdaE*p)*G.Gw)p,
                % Hprime = zMatrix*G.Gw + H*G.Guv
                LambdaE = tenzeros(3,2,2);
                LambdaE(1,1,1) = 1/2;  LambdaE(2,2,2) = 1/2; 
                LambdaE(3,1,2) = 1; 
                obj.initialization.LambdaE = LambdaE;
                                      
        end
            
     
        function [K,F] = tangent_stiffness_and_force(obj,p)
            
            %extract nodal displacements and nodal temeperatures from full
            %vector
            [pe] = obj.extract_element_data(p);
            
            %get the specifications for element integration
            X = obj.quadrature.X; %get quadrature points (natural coordinates)
            W = obj.quadrature.W; %get quadrature weights
            
            %initialize stiffness matrix and internal force vector
            K = obj.initialization.K; %tangent stiffness
            F = obj.initialization.F; %element internal force vector
            
            %material properites
            alpha = obj.initialization.alpha; %thermal expansion coefficient
            Cgamma = obj.initialization.Cgamma; %shear material constitutive law
            
            %retrieve strain displacements matricies
            H = obj.initialization.H;
            AEfun = obj.initialization.AEfun;
            AdelEfun = obj.initialization.AdelEfun;
            
            %retrieve membrane - bending constitutive matrices:
            Cm = obj.initialization.Cm;
            Cb = obj.initialization.Cb; 
            
            NintPoints = length(W);
            
            %start integration 
            for ii = 1:NintPoints
                
               xii = X(:,ii); %gauss point (natural coordinates)
               wii = W(ii); %integration weight
               
               %shape functions for the derivatives and detJ
               [shFun,G,detJ] = obj.Gfun(xii);
               
               dwdx = G.Gw*pe; %dwdx
               thxy = [shFun*obj.initialization.EMthx; %[thx;thy]
                         shFun*obj.initialization.EMthy];
               
               %compute the B matrices
               BE = H*G.Guv + AEfun(dwdx,G.dzdx)*G.Gw;
               Bgamma = thxy + G.Gw;
               Bchi = H*G.Gthxy;
               BdelE = H*G.Guv + AdelEfun(dwdx,G.dzdx)*G.Gw;
               
               %compute the internal actions
               N = Cm*(BE*pe); %membrane internal actions
               M = Cb*(Bchi*pe); %bending internal actions
               V = Cgamma*Bgamma*pe; %out of plane shear internal actions
               
               %integrate the internal force vector
               F = F + (BdelE'*N + Bchi'*M + Bgamma'*V)*detJ*wii;
               
               %compute the tangent stiffness matrix
               K1 = BdelE'*Cm*BdelE + Bchi'*Cb*Bchi + ...
                       + Bgamma'*Cgamma*Bgamma;
                
               K2 = (G.Gw).'*[N(1) N(3); 
                            N(3) N(2)]*G.Gw;
                
               K = K + (K1 + K2)*detJ*wii;  
                
            end          
            
        end
        
       
        
        function [T2, globalSubs] = T2(obj,varargin) 
            %tensor associated to polynomial forces of order 2
            
            if ~isempty(varargin)
                Ve = tensor(varargin{1});
                flag = true;
            else
                flag = false;
            end
            
            m = obj.nNodes*obj.nDOFPerNode; %number dofs per node
                
            % global DOFs associated to the element nodes
            index = get_index(obj.nodeIDs,obj.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = {index, index, index};                
              
            H = obj.initialization.H;
            % dE = (Hprime + (lambdaDelE*p)*G.Gw)dp, % E = (Hprime + (lambdaE*p)*G.Gw)p
            lambdaE = obj.initialization.LambdaE;
            lambdaDelE = obj.initialization.LambdaDelE;
            
            %retrieve membrane - bending constitutive matrices:
            Cm = obj.initialization.Cm;
            
            %get the specifications for element integration
            X = obj.quadrature.X; %get quadrature points (natural coordinates)
            W = obj.quadrature.W; %get quadrature weights
            NintPoints = length(W);
            
            %initialize output
            T2 = tenzeros(m,m,m); 
                          
            %start integration 
            for ii = 1:NintPoints
                
                xii = X(:,ii); %gauss point (natural coordinates)
                wii = W(ii); %integration weight

                %shape functions for the derivatives and detJ
                [~,G,detJ] = obj.Gfun(xii);
                
                Hprime = (H*G.Guv + obj.initialization.zMatrixFun(G.dzdx)*G.Gw);
                
                GwTen = tensor(G.Gw);
                LdelE = ttt(lambdaDelE,GwTen,3,1);
                LE = ttt(lambdaE,GwTen,3,1);
                
                dT21 = ttt( ttt(tensor(Hprime'*Cm),LE,2,1), GwTen,2,1) ;
                dT22 = permute( ttt( ttt(LdelE, GwTen, 2,1), tensor(Cm*Hprime), 1,1), [2,1,3]); 
                
                %integrate  
                T2 = T2 + (dT21 + dT22)*detJ*wii;    
                
            end
            
            if flag
                
                T2 = ttt(ttt(ttt(Ve,T2,1,1),Ve,2,1),Ve,2,1);
               
            end
                          
        end
        
              
        
        function [T3, globalSubs] = T3(obj,varargin) 
            %tensor associated to polynomial forces of order 3
            
            if ~isempty(varargin)
                Ve = tensor(varargin{1});
                flag = true;
            else
                flag = false;
            end
            
            m = obj.nNodes*obj.nDOFPerNode; %number dofs per node
                
            % global DOFs associated to the element nodes
            index = get_index(obj.nodeIDs,obj.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = {index, index, index, index}; 
                
            % dE = (Hprime + (lambdaDelE*p)*G.Gw)dp, 
            lambdaE = obj.initialization.LambdaE;
            
            % E = (Hprime + (lambdaE*p)*G.Gw)p
            lambdaDelE = obj.initialization.LambdaDelE;
            
            %retrieve membrane - bending constitutive matrices:
            Cm = obj.initialization.Cm;
            
            %get the specifications for element integration
            X = obj.quadrature.X; %get quadrature points (natural coordinates)
            W = obj.quadrature.W; %get quadrature weights
            NintPoints = length(W);
            
            %initialize output
            T3 = tenzeros(m,m,m,m); 
                          
            %start integration 
            for ii = 1:NintPoints
                
                xii = X(:,ii); %gauss point (natural coordinates)
                wii = W(ii); %integration weight

                %shape functions for the derivatives and detJ
                [~,G,detJ] = obj.Gfun(xii);
                
                GwTen = tensor(G.Gw);
                LdelE = ttt(lambdaDelE,GwTen,3,1);
                LE = ttt(lambdaE,GwTen,3,1);
                
                T31 = ttt(LdelE,GwTen,2,1);
                T32 = ttt(LE,GwTen,2,1);
                T33 = ttt(tensor(Cm),T32,2,1);
                
                dT3ii = permute( ttt(T31,T33,1,1), [2,1,3,4]);
                
                T3 = T3 + dT3ii*detJ*wii; %integrate     
                
            end         
            
            if flag
                
                T3 = ttt(ttt(ttt(ttt(Ve,T3,1,1),Ve,2,1),Ve,2,1),Ve,2,1);
               
            end
            
        end
            
                    

        function Kd = stiffness_derivative(obj, v, varargin)
            
            % this function compute the derivative of the tangent stiffness 
            % constructed around point p, along the direction v
            
            %p is the linearization point
            %v is the expansion direction        

            index = get_index(obj.nodeIDs,obj.nDOFPerNode); 
            
            if isempty(varargin)
                pe = zeros(obj.nNodes*obj.nDOFPerNode,1);
            else
                p = varargin{1};
                pe = p(index,:); % get the nodal displacements from global displacements
            end
            
            ve = v(index,:); % get the nodal directions from global direction
            
            %get the specifications for element integration
            X = obj.quadrature.X; %get quadrature points (natural coordinates)
            W = obj.quadrature.W; %get quadrature weights
            
            %retrieve membrane - bending constitutive matrices:
            Cm = obj.initialization.Cm;
            
            %retrieve strain displacements matricies
            H = obj.initialization.H;
            AdelEfun = obj.initialization.AdelEfun;
            
            % dE = (Hprime + (lambdaDelE*p)*G.Gw)dp, 
            % E = (Hprime + (lambdaE*p)*G.Gw)p
            lambdaDelE = obj.initialization.LambdaDelE;
            
            %initialize stiffness derivative matrix
            Kd = obj.initialization.K; %tangent stiffness
            
            NintPoints = length(W);
            
            %start integration 
            for ii = 1:NintPoints
                
               xii = X(:,ii); %gauss point (natural coordinates)
               wii = W(ii); %integration weight
               
               %shape functions for the derivatives and detJ
               [~,G,detJ] = obj.Gfun(xii);
               
               dwdx = G.Gw*pe; %dwdx
               
               %compute the B matrix
               BdelE = H*G.Guv + AdelEfun(dwdx,G.dzdx)*G.Gw;
               
               %compute delta of membrane strain
               GwTen = tensor(G.Gw);
               LdelE = ttt(lambdaDelE,GwTen,3,1);
               
               dBdelE = double(ttv(LdelE,ve,3))*G.Gw;
               
               %dK1
               dK1ii =  dBdelE'*Cm*BdelE + BdelE'*Cm*dBdelE;
               
               %compute dN
               dN = Cm*BdelE*ve; 
               
               %dK2
               dK2ii = G.Gw'*[dN(1) dN(3); dN(3) dN(2)]*G.Gw;
               
               %integrate
               Kd = Kd + (dK1ii + dK2ii)*detJ*wii;             
                
            end
               
        end
        
        
                        
        function [xe]  = extract_element_data(obj,x)    
            % this function extracts the nodal displacements and T from 
            % full vectors
            
            % x is a vector of full DOFs
            % T is the vector of nodal temperatures
            index = get_index(obj.nodeIDs,obj.nDOFPerNode);
            
            xe = x(index,:); % get the element displacements from global displacements

        end
        
        
        
        function M = mass_matrix(obj)
            
            t = obj.thickness;
            pho = obj.Material.DENSITY; 
            
            EMu = obj.initialization.EMu;
            EMv = obj.initialization.EMv;
            EMw = obj.initialization.EMw;
            EMthx = obj.initialization.EMthx;
            EMthy = obj.initialization.EMthy;
            
            %get the specifications for element integration
            X = obj.quadrature.X; %get quadrature points (natural coordinates)
            W = obj.quadrature.W; %get quadrature weights
            
            NintPoints = length(W); %number of integration points
            
            zeri = zeros(1,size(EMu,2));
            
            M = obj.initialization.M; %mass matrix
            
            for ii = 1:NintPoints
                
            xii = X(:,ii); %gauss point (natural coordinates)
            wii = W(ii); %integration weight
            [N,detJ,~] = shape_function_derivatives(obj,xii);
            
            huvw = [N*EMu; N*EMv; N*EMw];
            hthxy = [N*EMthx; N*EMthy; zeri];
            
            M = M + (huvw'*huvw*t + hthxy'*hthxy*t^3/12)*pho*detJ*wii;
            
            end
            
        end

        function Nfield = membrane_forces(obj,nodeIDsQueried,p)

            %extract element data from full vector of unknowns
            pe  = extract_element_data(obj,p);

            %extract membrane forces at gauss points
            %get the specifications for element integration
            X = obj.quadrature.X; %get quadrature points (natural coordinates)
            
            %material properites
            alpha = obj.initialization.alpha; %thermal expansion coefficient
            Cgamma = obj.initialization.Cgamma; %shear material constitutive law
            
            %retrieve strain displacements matricies
            H = obj.initialization.H;
            AEfun = obj.initialization.AEfun;
            
            %retrieve membrane - bending constitutive matrices:
            Cm = obj.initialization.Cm;
            Cb = obj.initialization.Cb; 
            
            NintPoints = size(X,2);
            
            %initialize membrane force vector at Gauss Points (row -> gauss
            %point, column -> Nx, Ny, Nxy
            membraneForcesGauss = zeros(NintPoints,3); %[Nx, Ny, Nxy]
            

            for ii = 1:NintPoints

               xii = X(:,ii); %gauss point (natural coordinates)
               
               %shape functions for the derivatives and detJ
               [shFun,G,detJ] = obj.Gfun(xii); %if you initialize this object, the memory req. increase, but you are faster online
               
               dwdx = G.Gw*pe; %dwdx
               thxy = [shFun*obj.initialization.EMthx; %[thx;thy]
                         shFun*obj.initialization.EMthy];
               
               %compute the B matrices
               BE = H*G.Guv + AEfun(dwdx,G.dzdx)*G.Gw;
               Bgamma = thxy + G.Gw;
               Bchi = H*G.Gthxy;
               
               %compute the internal actions
               membraneForcesGauss(ii,:) = (Cm*BE*pe)'; %membrane internal actions

            end
            
            %extrapolate value of field from gauss points to nodes
            gauss2NodesMatrix = obj.gauss_to_nodes();
            gauss2QueriedNodesMatrix = gauss2NodesMatrix(ismember(obj.nodeIDs,nodeIDsQueried),:);
            Nfield = gauss2QueriedNodesMatrix*membraneForcesGauss;
            
        end


        function Mbfield = bending_moments(obj,nodeIDsQueried,p)

            %extract element data from full vector of unknowns
            pe  = extract_element_data(obj,p);

            %extract membrane forces at gauss points
            %get the specifications for element integration
            X = obj.quadrature.X; %get quadrature points (natural coordinates)
            
            %material properites
            alpha = obj.initialization.alpha; %thermal expansion coefficient
            
            %retrieve strain displacements matricies
            H = obj.initialization.H;
            
            %retrieve bending constitutive matrices:
            Cb = obj.initialization.Cb; 
            
            NintPoints = size(X,2);
            
            %initialize membrane force vector at Gauss Points (row -> gauss
            %point, column -> Mx, My, Mxy
            momentsGauss = zeros(NintPoints,3); %[Mx, My, Mxy]
            

            for ii = 1:NintPoints

               xii = X(:,ii); %gauss point (natural coordinates)
               
               %shape functions for the derivatives and detJ
               [shFun,G,detJ] = obj.Gfun(xii); %if you initialize this object, the memory req. increase, but you are faster online
               
               %compute the B matrices
               Bchi = H*G.Gthxy;
               
               %compute the internal actions
               momentsGauss(ii,:) = (Cb*Bchi*pe)'; %membrane internal actions

            end
            
            %extrapolate value of field from gauss points to nodes
            gauss2NodesMatrix = obj.gauss_to_nodes();
            gauss2QueriedNodesMatrix = gauss2NodesMatrix(ismember(obj.nodeIDs,nodeIDsQueried),:);
            Mbfield = gauss2QueriedNodesMatrix*momentsGauss;
            
        end



        function Vfield = shear_forces(obj,nodeIDsQueried,p)

            %extract element data from full vector of unknowns
            pe  = extract_element_data(obj,p);

            %extract membrane forces at gauss points
            %get the specifications for element integration
            X = obj.quadrature.X; %get quadrature points (natural coordinates)
            
            %material properites
            Cgamma = obj.initialization.Cgamma; %shear material constitutive law
 
            NintPoints = size(X,2);
            
            %initialize shear force vector at Gauss Points (row -> gauss
            %point, column -> Vx, Vy
            shearForcesGauss = zeros(NintPoints,2); %[Vx, Vy] %shear force (z direction)
            

            for ii = 1:NintPoints

               xii = X(:,ii); %gauss point (natural coordinates)
               
               %shape functions for the derivatives and detJ
               [shFun,G,detJ] = obj.Gfun(xii); %if you initialize this object, the memory req. increase, but you are faster online
               
               thxy = [shFun*obj.initialization.EMthx; %[thx;thy]
                         shFun*obj.initialization.EMthy];
               
               %compute the B matrix
               Bgamma = thxy + G.Gw;
               
               %compute out of plane shear internal actions
               shearForcesGauss(ii,:) = (Cgamma*Bgamma*pe)'; %membrane internal actions

            end
            
            %extrapolate value of field from gauss points to nodes
            gauss2NodesMatrix = obj.gauss_to_nodes();
            gauss2QueriedNodesMatrix = gauss2NodesMatrix(ismember(obj.nodeIDs,nodeIDsQueried),:);
            Vfield = gauss2QueriedNodesMatrix*shearForcesGauss;
            
        end
        
        
           
        function setZ(obj,znodes) 
          %function to store the values of z in object
          %znodes contains the z height of the nodes (full vector)
            obj.z = znodes(obj.nodeIDs);
            
        end
        
        
       
        function [N,Gout,detJ] = Gfun(obj,X)
            
        %this function computes the derivatives from nodal displacements and 
        %shape functions through matrix G
        
            %X -> vector of natural coordinates
            %p -> vector of nodal displacements (in the local ref. frame)
            
            [N,detJ,dNdx] = shape_function_derivatives(obj,X); 
            
            %in plane disp,  Guv*p -> [dudx,dudy,dvdx,dvdy]'
            Gout.Guv = [dNdx * (obj.initialization.EMu);
                dNdx * (obj.initialization.EMv)]; 
            
            %out of plane disp,   Gw*p -> [dwdx,dwdy]'
            Gout.Gw = (dNdx * (obj.initialization.EMw)); 
            
            %bending rotations,  Gthxy*p -> [dthxdx,dthxdy,dthydx,dthydy]'
            Gout.Gthxy = [dNdx * (obj.initialization.EMthx);
                dNdx * (obj.initialization.EMthy)]; 
            
            %z coordinates
            Gout.dzdx = dNdx * obj.z;  % dzdx = [dzdx;dzdy]   
            
        end
        
              
       function quadrature_rule(obj, Ngauss)
        %this function defines the quadrature rule for the element
                DIM = obj.nDim;
                ELTYPE = obj.elType;
                switch ELTYPE
                    case {'QUAD'}
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
                        obj.quadrature.Ng = Ngauss;
                        obj.quadrature.X = X;	% gauss integration points
                        obj.quadrature.W = W;	% gauss integration weights
%                     case {'TET'} %CHECK IF IT IS SUITED FOR TRIANGULAR
%                     ELEMENTS!
%                         [x, w] = inttet(Ngauss);
%                         self.quadrature.Ng = Ngauss;
%                         self.quadrature.X = x;	% gauss integration points
%                         self.quadrature.W = w;	% gauss integration weights
%                     case {'WED'}
%                         [w, x] = wedge_rule(Ngauss.lin, Ngauss.tri);
%                         self.quadrature.Ng = Ngauss;
%                         self.quadrature.X = x;	% gauss integration points
%                         self.quadrature.W = w;	% gauss integration weights
                    otherwise
                        error([' No quadrature rule found for the element "' ...
                            ELTYPE '". Accepted values are: HEX/QUAD, TET/TRI, WED'])
                end
       end
                
       
       
    end
    
    
    
end

