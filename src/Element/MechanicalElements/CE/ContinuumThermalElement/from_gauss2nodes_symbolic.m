
sr3 = sqrt(3);

X = sr3*[-1    -1    -1;
         1    -1    -1;
         1     1    -1;
        -1     1    -1;
        -1    -1     1;
         1    -1     1;
         1     1     1;
        -1     1     1;
         0    -1    -1;
         1     0    -1;
         0     1    -1;
        -1     0    -1;
         0    -1     1;
         1     0     1;
         0     1     1,
        -1     0     1;
        -1    -1     0;
         1    -1     0;
         1     1     0;
        -1     1     0];

nNodes = 20;
N = zeros(nNodes,8);

for ii = 1:nNodes
    xi = X(ii,1);
    eta = X(ii,2);
    zeta = X(ii,3);
    N(ii,:) = shape_functions(xi,eta,zeta);
end


function N = shape_functions(g,h,r)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 8-NODED HEXAHEDRON
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.4 Solid isoparametric quadrilaterals and hexahedra
          
            Nor = 1/8*[(1-g)*(1-h)*(1-r); 
                (1+g)*(1-h)*(1-r);
                     (1+g)*(1+h)*(1-r); 
                     (1-g)*(1+h)*(1-r);
                     (1-g)*(1-h)*(1+r); 
                     (1+g)*(1-h)*(1+r);
                     (1+g)*(1+h)*(1+r);
                     (1-g)*(1+h)*(1+r)];

            N = zeros(1,8);
            N(1:end) = [Nor(7),Nor(3),Nor(6),Nor(2),Nor(8),Nor(4),Nor(5),Nor(1)];
end

