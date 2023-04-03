function err = residual_error_ROM(residualFun,ROMdisp,ROMvel,ROMacc,ROMtime)

%this function computes the norm of the residual of the ROM solution 
%over time

nTimeSamples = length(ROMtime);

err = zeros(nTimeSamples,1);

for ii = 1:nTimeSamples
  
  t = ROMtime(ii);
  
  q = ROMdisp(:,ii);
  qd = ROMvel(:,ii);
  qdd = ROMacc(:,ii);
  
  err(ii) = norm(residualFun(q,qd,qdd,t));

end


end