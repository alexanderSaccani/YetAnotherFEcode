function Vt = V_fun_t(array4int,psamples,pfunt,tmin,tmax,ntsamples)

time = linspace(tmin,tmax,ntsamples);

%nROMs = size(array4int,1);
nVecBasis = size(array4int,3);
nDofsFull = size(array4int,2);

Vt = zeros(nDofsFull,nVecBasis);

xsamples = psamples(1,:);
for jj = 1:ntsamples

    pCurr = pfunt(time(jj));
    xCurr = pCurr(1); 

    parfor ii = 1:nVecBasis
    
    Vt(:,ii,jj) = interp1(xsamples,squeeze(array4int(:,:,ii)),xCurr,'spline');
    
    end

end






end