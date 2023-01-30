function disp = quasistatic_solution( Assembly, T, time )

%T -> temperature function handle of time t: @(t)T 
%time -> vector of queried time samples

n_time_samples = length(time); %number of time samples
nDofsF = Assembly.Mesh.EBC.nDOFs; % #dofs free (without bc)

%forcing = 0
F = zeros(nDofsF,1);

%initialize guess
u_pre = zeros(nDofsF,1);

disp = zeros(nDofsF,n_time_samples); %initialize output
for ii = 1:n_time_samples
    t = time(ii);
    [~,u_ii] = static_equilibrium_thermal(Assembly, F, T(t),'initialGuess', u_pre);
    disp(:,ii) = u_ii; %save output
    u_pre = u_ii; %update initial guess
end

end