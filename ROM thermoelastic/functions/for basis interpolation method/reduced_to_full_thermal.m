function u_full = reduced_to_full_thermal(red_solution,time,ROMs,p_th,u_quasist,time_quasist)
%REDUCED_TO_FULL_THERMAL Summary of this function goes here
%   Detailed explanation goes here

n_t_samples = length(time);

q = red_solution;

%initialize output
u_full = zeros(size(ROMs.models{1}.V,1), n_t_samples);

% full assembly (for interpolation with quasistatic solution)
fullAssembly = ROMs.fullAssembly;

for ii = 1:n_t_samples
    
    t = time(ii);
    
    %interpolate basis according to time istant
    %sample the thermal parameter
    p = p_th(t);

    %compare the thermal parameter with sampled parameters (for which ROM is
    %available)
    p_th_samples = ROMs.parameters; %vector of sampled parameters for which models in ROMs were constructed

    [~,ind_p1] = min(abs(p_th_samples - p));
    p_th_samples_tmp = p_th_samples;
    p_th_samples_tmp(ind_p1) = inf;
    [~,ind_p2] = min(abs(p_th_samples_tmp - p));

    p1 = p_th_samples(ind_p1); %p, thermal parameter follows in the interval of bounds p1 and p2 (not necesserely in order)
    p2 = p_th_samples(ind_p2);

    %interpolate basis
    V1 = ROMs.models{ind_p1}.V;
    V2 = ROMs.models{ind_p2}.V;
    V = V1+(p-p1)/(p2-p1)*(V2-V1); %linear interpolation of basis
% 
%     %%interpolation of thermal equilibrium with spline from samples of constructed ROMs does not work
%     u_eq1 = ROMs.models{ind_p1}.thermal_eq;
%     u_eq2 = ROMs.models{ind_p2}.thermal_eq;
%     u_eq = u_eq1+(p-p1)/(p2-p1)*(u_eq2-u_eq1);
%     
%     %interpolation of thermal equilibrium with spline from samples of constructed ROMs does not work
%     eqq = zeros(size(V,1),length(ROMs.parameters));
%     for jj = 1:length(ROMs)
%         eqq(:,jj) =  ROMs.models{jj}.thermal_eq;
%     end
%     u_eq = spline(p_th_samples,eqq,p);
% 

    % interpolate the thermal equilibrium from a quasistatic full solution
    u_eq = interp1(time_quasist.',u_quasist.',t);
    u_eq = u_eq.';
    u_eq = fullAssembly.constrain_vector(u_eq);


    u_full(:,ii) = u_eq + V*q(:,ii);


end

end
