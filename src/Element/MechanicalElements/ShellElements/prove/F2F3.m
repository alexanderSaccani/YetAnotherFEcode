%         function F3 = F3(obj,p,T,gradT)
%             
%             [pe,Te,gradTe] = obj.extract_element_data(p,T,gradT);
%             
%             F3 = zeros(length(pe),1);
%             Cm = obj.initialization.Cm;
%             Cb = obj.initialization.Cb; 
%             
%             %get the specifications for element integration
%             X = obj.quadrature.X; %get quadrature points (natural coordinates)
%             W = obj.quadrature.W; %get quadrature weights
%             NintPoints = length(W);
%             
%             lambdaE = obj.initialization.LambdaE;
%             lambdaDelE = obj.initialization.LambdaDelE;
%             
%             for ii = 1:NintPoints
%                 
%                 xii = X(:,ii); %gauss point (natural coordinates)
%                 wii = W(ii); %integration weight
% 
%                 %shape functions for the derivatives and detJ
%                 [shFun,G,detJ] = obj.Gfun(xii);
%                 
% %               Hprime = (H*G.Guv + obj.initialization.zMatrixFun(G.dzdx)*G.Gw);
%                 GwTen = tensor(G.Gw);
%                 LdelE = ttt(lambdaDelE,GwTen,3,1);
%                 LE = ttt(lambdaE,GwTen,3,1);
%                 
%                 df3ii = (double(ttv(LdelE,pe,3))*G.Gw).'*Cm*double(ttv(LE,pe,3))*G.Gw*pe;
%                 F3 = F3 + df3ii*detJ*wii;
%    
%                 
%             end
%             
%             
%         end
%         
%         
%         function F2 = F2(obj,p,T,gradT)
%             
%             [pe,Te,gradTe] = obj.extract_element_data(p,T,gradT);
%             
% 
%             Cm = obj.initialization.Cm;
%             Cb = obj.initialization.Cb; 
%             H = obj.initialization.H;
%             
%             
%             %get the specifications for element integration
%             X = obj.quadrature.X; %get quadrature points (natural coordinates)
%             W = obj.quadrature.W; %get quadrature weights
%             NintPoints = length(W);
%             
%             lambdaE = obj.initialization.LambdaE;
%             lambdaDelE = obj.initialization.LambdaDelE;
%             
%             
%             F2 = zeros(length(pe),1);
%             for ii = 1:NintPoints
%                 
%                 xii = X(:,ii); %gauss point (natural coordinates)
%                 wii = W(ii); %integration weight
% 
%                 %shape functions for the derivatives and detJ
%                 [shFun,G,detJ] = obj.Gfun(xii);
%                 
%                 Hprime = (H*G.Guv + obj.initialization.zMatrixFun(G.dzdx)*G.Gw);
%                 GwTen = tensor(G.Gw);
%                 LdelE = ttt(lambdaDelE,GwTen,3,1);
%                 LE = ttt(lambdaE,GwTen,3,1);
%                 
%                 df2ii = (Hprime'*Cm)*(double(ttv(LE,pe,3))*G.Gw)*pe + ...
%                     ( double(ttv(LdelE,pe,3)) *G.Gw)'*(Cm*Hprime)*pe;
%                 
%                 F2 = F2 + df2ii*detJ*wii;
%    
%                 
%             end
%             
%             
%         end
