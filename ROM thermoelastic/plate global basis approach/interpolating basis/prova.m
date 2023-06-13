% prova interpolazione
clear all
close all
clc

%%
tmin = 0;
tmax = 1;

%psamples = [1,3,4,5,6.7,9,10];
psamples = linspace(1,10,10);
psamples(2) = 1;
psamples(1) = 2;
ysamples =  [2,4,5,10,14,10,9,12,11,13];


t = linspace(0,tmax,100);
interpolatedy = zeros(1,length(t));
pt = @(t) pfunt(t);

for ii = 1:length(t)
    tt = t(ii);
    interpolatedy(ii) = interpolate(psamples',ysamples,pt(tt));
    
end


figure
plot(psamples,ysamples,'*')
hold on
plot(pt(t),interpolatedy)
%%
function out = pfunt(t)
 out = 10*t;
end

function out = interpolate(psamples,ysamples,pcurr)
 

    nSamples = size(psamples,2);
    dist_p_psamples1 = abs(psamples - repmat(pcurr,1,nSamples));
    [dist1,closestModel1] = min(dist_p_psamples1);
    y1 = ysamples(closestModel1);
    p1 = psamples(closestModel1);
    x1 = p1;
    
    disp(['closest model is: ',num2str(closestModel1)])

    psamples(closestModel1) = psamples(closestModel1)*inf;
    dist_p_psamples2 = abs(psamples - repmat(pcurr,1,nSamples));
    [dist2,closestModel2] = min(dist_p_psamples2);
    y2 = ysamples(closestModel2);
    p2 = psamples(closestModel2);
    x2 = p2;

    disp(['closest model is: ',num2str(closestModel2)])
    
    %interpolate the basis
    %V = dist2/distp1p2*V1 + dist1/distp1p2*V2;
    out = y1 + (pcurr-x1)/(x2-x1)*(y2-y1);


end
