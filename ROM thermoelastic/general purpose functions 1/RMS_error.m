function err = RMS_error(HFMsol,HFMtime,ROMsol,ROMtime)

% Author Alexander Saccani, phd candidate ETH Zurich, 4/2023

%this function computes the RMS error of the ROM solution w.r.t. the HFM
%solution.
%
% INPUT: HFMsol: double array in which the columns correspond to different
%               time samples of HFM sol
%        HFMtime: 1D double array in which are stored the time samples of
%               the HFM analysis
%        ROMsol: same as HFMsol
%        ROMtime: same as HFMtime
% OUTPUT: time history of RMS, referred to HFM time samples

ROMsol = (interp1(ROMtime',ROMsol',HFMtime'))';

error = HFMsol - ROMsol;
sqrdErr = error.^2;

meanSol = sum(abs(HFMsol),1)/size(HFMsol,1); %normalized with respect to mean solution
err = sqrt(sum(sqrdErr,1))/(size(error,1))./meanSol;

end








