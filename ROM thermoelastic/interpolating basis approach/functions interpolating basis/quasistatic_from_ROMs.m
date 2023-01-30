function u_qs = quasistatic_from_ROMs(time,ROMs,p_th)


    n_time_sampl = length(time);
    n_dof = size(ROMs.models{1}.thermal_eq,1);
    
    p_th_samples = ROMs.parameters; %vector of sampled parameters for which models in ROMs were constructed
    
    u_qs = zeros(n_dof,n_time_sampl);

    for ii = 1:n_time_sampl
        
        p = p_th(time(ii));

        [~,ind_p1] = min(abs(p_th_samples - p));
        p_th_samples_tmp = p_th_samples;
        p_th_samples_tmp(ind_p1) = inf;
        [~,ind_p2] = min(abs(p_th_samples_tmp - p));

        p1 = p_th_samples(ind_p1); %p, thermal parameter follows in the interval of bounds p1 and p2 (not necesserely in order)
        p2 = p_th_samples(ind_p2);

        %interpolate full displacement vector at equilibrium
        u_eq1 = ROMs.models{ind_p1}.thermal_eq;
        u_eq2 = ROMs.models{ind_p2}.thermal_eq;
        u_qs(:,ii) = u_eq1+(p-p1)/(p2-p1)*(u_eq2-u_eq1);
        
    end

end