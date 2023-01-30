function [ r, drdqdd,drdqd,drdq, c0] = residual_ROM_thermal( q, qd, qdd, t, ROMs, Fext, p_th, T_fn_p, u_quasist, time_quasist)

% INPUT: p_th is a function handle of time that s
%        T_fn_p is a function handle of p_thermal, that describes how the
%        nodal temperatures vary with p_thermal
%        ROMs should have fields models, parameters, and fullAssembly
%        u_quasistatic, is the full quasistatic solution of the dynamic
%        problem
%        time_quasist is the time vector of the provided quasistatic solution 

% full assembly
fullAssembly = ROMs.fullAssembly;

%% Interpolate the basis and the equilibrium point
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

%compute T(p_th(t))
T = T_fn_p(p);

%interpolate basis
V1 = ROMs.models{ind_p1}.V;
V2 = ROMs.models{ind_p2}.V;
V = V1+(p-p1)/(p2-p1)*(V2-V1); %linear interpolation of basis

% Interpolate the thermal equilibrium position from constructed ROMs
% doesn't work
% interpolation with spline
% eqq = zeros(size(V,1),length(ROMs.parameters));
% for ii = 1:length(ROMs)
%     eqq(:,ii) =  ROMs.models{ii}.thermal_eq;
% end
% u_eq = spline(p_th_samples,eqq,p);

% %interpolate full displacement vector at equilibrium
% %linear interpolation from sampled ROMs does not work
% u_eq1 = ROMs.models{ind_p1}.thermal_eq;
% u_eq2 = ROMs.models{ind_p2}.thermal_eq;
% u_eq = u_eq1+(p-p1)/(p2-p1)*(u_eq2-u_eq1);

% interpolate the thermal equilibrium from a quasistatic full solution
u_eq = interp1(time_quasist.',u_quasist.',t);
u_eq = u_eq.';
u_eq = fullAssembly.constrain_vector(u_eq);

%% compute ROM
M = fullAssembly.constrain_matrix(fullAssembly.DATA.M); %can be more efficient by providing the constrain matrix (then pay attention to ROMs)
C = fullAssembly.constrain_matrix(fullAssembly.DATA.C);

Mr = V.'*M*V; %constrain matrix
Cr = V.'*C*V; %constrain matrix

u = fullAssembly.unconstrain_vector(u_eq + V*q); %get full coordinates

[Ktg, Fi_full] = fullAssembly.tangent_stiffness_and_force(u,T); %get full internal force and tangent stiffness
F_int_r = V.'*fullAssembly.constrain_vector(Fi_full); %back to reduced space
Ktgr = V.'*fullAssembly.constrain_matrix(Ktg)*V;

%compute residual
F_inertial = Mr*qdd;
F_damping = Cr*qd;
F_ext_r = V.'*fullAssembly.constrain_vector(Fext(t));

r = F_inertial + F_damping + F_int_r - F_ext_r;

%Jacobians
drdqdd = Mr;
drdqd = Cr;
drdq = Ktgr;
%% 
% We use the following measure to comapre the norm of the residual $\mathbf{r}$
% 
% $$\texttt{c0} = \|\mathbf{M_V}\ddot{\mathbf{q}}\| + \|\mathbf{C_V}\dot{\mathbf{q}}\| 
% + \|\mathbf{F_V}(\mathbf{q})\| + \|\mathbf{V}^{\top}\mathbf{F}_{ext}(t)\|$$
c0 = norm(F_inertial) + norm(F_damping) + norm(F_int_r) + norm(F_ext_r);
end