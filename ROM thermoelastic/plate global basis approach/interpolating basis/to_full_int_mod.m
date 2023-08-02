function fullDisp = to_full_int_mod(dispRed,timeRed,pfunt,ROMs,array4int)

numb_time_samples = size(dispRed,2);
ndofsFull = size(ROMs.models{1}.V,1);

fullDisp = zeros(ndofsFull,numb_time_samples);
psamples = ROMs.psamples;
xsamples = psamples(1,:);

parfor ii = 1:numb_time_samples
    
    t = timeRed(ii);
    %scegli la base

    %interpolate the basis
    p_current = pfunt(t);
    xc = p_current(1);
    
   
% 
%     %p is column vector
%     nSamples = size(psamples,2);
%     dist_p_psamples1 = sqrt(sum((psamples - repmat(p_current,1,nSamples)).^2));
%     [dist1,closestModel1] = min(dist_p_psamples1);
%     V1 = ROMs.models{closestModel1}.V;
%     p1 = psamples(:,closestModel1);
%     x1 = p1(1);
%     
%     disp(['model used is: ',num2str(closestModel1)])
%     
%     psamples(:,closestModel1) = psamples(:,closestModel1)*inf;
%     dist_p_psamples2 = sqrt(sum((psamples - repmat(p_current,1,nSamples)).^2));
%     [dist2,closestModel2] = min(dist_p_psamples2);
%     V2 = ROMs.models{closestModel2}.V;
%     p2 = psamples(:,closestModel2);
%     x2 = p2(1);
%     
%     disp(['model used is: ',num2str(closestModel2)])
%     
%     distp1p2 = sqrt(sum((p1-p2).^2));
%     
%     %interpolate the basis
%     %V = dist2/distp1p2*V1 + dist1/distp1p2*V2;
%     V = V1 + (xc-x1)/(x2-x1)*(V2-V1);
%% interpolate the basis

nROMs = size(array4int,1);
nVecBasis = size(array4int,3);
nDofsFull = size(array4int,2);

V = zeros(nDofsFull,nVecBasis);


for jj = 1:nVecBasis

a = squeeze(array4int(:,:,jj));
V(:,jj) = interp1(xsamples,a,xc,'spline');

end

%store the full solution
fullDisp(:,ii) = V*dispRed(:,ii);

end


end