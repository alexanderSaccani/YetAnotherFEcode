%anticlockwise numeration of nodes, first external nodes and then midpoint
%nodes
sr3 = sqrt(3);

nNodes = 4;

if nNodes == 4
    X = [-sr3,-sr3;
        +sr3,-sr3;
        sr3,sr3;
        -sr3,sr3];
elseif nNodes == 8
    X = [-sr3,-sr3;
        +sr3,-sr3;
        sr3,sr3;
        -sr3,sr3;
        0,-sr3;
        sr3,0;
        0,sr3;
        -sr3,0];
end

N = zeros(nNodes,4);

for ii = 1:nNodes
    xi = X(ii,1);
    eta = X(ii,2);
    N(ii,:) = shape_fun(xi,eta);
end


% original shape fun
% function fun = shape_fun(xi,eta)
%              fun = 1/4*[ +(xi - 1)*(eta - 1); %shape function (1d)
%                       -(xi + 1)*(eta - 1);
%                       +(xi + 1)*(eta + 1); 
%                       -(xi - 1)*(eta + 1)].';   
% end

function fun = shape_fun(xi,eta)
%for gauss points ([(1/sr3,1/sr3), (1/sr3,-1/sr3), (-1/sr3,1/sr3), (-1/sr3,-1/sr3)])
             fun = 1/4*[ +(xi + 1)*(eta + 1); %shape function (1d)
                      -(xi + 1)*(eta - 1);
                      -(xi - 1)*(eta + 1); 
                      +(xi - 1)*(eta - 1)].';   
end

% function fun = shape_fun(xi,eta)
%              fun = 1/4*[ -(xi - 1)*(eta + 1); %shape function (1d)
%                       -(xi + 1)*(eta - 1);
%                        +(xi - 1)*(eta - 1); 
%                       +(xi + 1)*(eta + 1)].';   
% end


% function fun = shape_fun(xi,eta)
%         fun = 1/4*[...
%                 -(1-xi)*(1-eta)*(1+xi+eta);  %shape fun 1D
%                 -(1+xi)*(1-eta)*(1-xi+eta);
%                 -(1+xi)*(1+eta)*(1-xi-eta); 
%                 -(1-xi)*(1+eta)*(1+xi-eta);
%                 2*(1-xi)*(1+xi)*(1-eta);  
%                 2*(1-eta)*(1+eta)*(1+xi);
%                 2*(1-xi)*(1+xi)*(1+eta);  
%                 2*(1-eta)*(1+eta)*(1-xi)].';  
% end