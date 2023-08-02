function quasistSol = interpolated_quasistatic(timeRedInt,psamples,equilibria,pfunt)

TimeSamples = 1:1:length(timeRedInt);

neq = length(equilibria);
nDofsFull = size(equilibria{1},1);
array4Int = zeros(neq,nDofsFull);
parfor ii = 1:neq

    array4Int(ii,:) = equilibria{ii}';

end

xsamples = psamples(1,:);

quasistSol = zeros(nDofsFull,length(timeRedInt));
parfor ii = TimeSamples

   tii = timeRedInt(ii);
   pii = pfunt(tii);

   quasistSol(:,ii) =  interp1(xsamples,array4Int,pii(1),'spline');

end

end