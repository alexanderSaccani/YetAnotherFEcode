function err = RMS_error_ROM(HFMdisp,HFMtime,ROMdisp,ROMtime,gdl,nodes)

% Author Alexander Saccani, phd candidate ETH Zurich, 4/2023

%this function computes the RMS error of the ROM solution w.r.t. the HFM
%solution.
%
% INPUT: HFMdisp: 3D double array in which the indeces correspond to (node, dof of node, tsamp)
%        HFMtime: 1D double array in which are stored the time samples of
%        the HFM analysis
%        ROMdisp: same as HFMdisp
%        ROMtime: same as HFMtime
%        nodes: indeces of the nodes for which that contribute to the
%        error. 
%        gdl: vector containing the number of dofs over which the error is
%        computed.
%
% OUTPUT: time history of RMS, referred to HFM time samples

HFMofInt = HFMdisp(nodes,gdl,1:end);
ROMofInt = ROMdisp(nodes,gdl,1:end);

nTimeSamplesHFM = length(HFMtime);
nTimeSamplesROM = length(ROMtime);

nNodes = size(HFMofInt,1);
nGdl = size(HFMofInt,2);

HFMsamples = zeros(nNodes*nGdl,nTimeSamplesHFM);
ROMsamples = zeros(nNodes*nGdl,nTimeSamplesROM);

for ii = 1:nTimeSamplesHFM
   
    HFMsamples(:,ii) = reshape(squeeze(HFMofInt(:,:,ii)),[],1); 
    
end

for ii = 1:nTimeSamplesROM
   
    ROMsamples(:,ii) = reshape(squeeze(ROMofInt(:,:,ii)),[],1);
    
end

ROMsamples = (interp1(ROMtime',ROMsamples',HFMtime'))';

error = HFMsamples - ROMsamples;
sqrdErr = error.^2;

%meanDisp = sum(abs(HFMsamples),1)/size(HFMsamples,1);
meanDisp = sum((max(HFMsamples,[],2) - min(HFMsamples,[],2)),1)/(nNodes*nGdl); %take into account time history of displacements to compute the normalization factor
err = sqrt(sum(sqrdErr,1))/(size(error,1))./meanDisp;



end








