function [ r, drdqdd,drdqd,drdq, c0] = residual_ROM_thermal_constBasis( q, qd, qdd, t, ROMs, Fext, p_th, T_fn_p)

% INPUT: p_th is a function handle of time that s
%        T_fn_p is a function handle of p_thermal, that describes how the
%        nodal temperatures vary with p_thermal
%        ROMs should have fields models, parameters, and fullAssembly
%        N.B.: in this functions in the first models, it must be stored the
%        basis corresponding to the cold structure. 

%% Interpolate the thermal equilibrium point
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

% %interpolate basis
% V1 = ROMs.models{ind_p1}.V;
% V2 = ROMs.models{ind_p2}.V;
% V = V1+(p-p1)/(p2-p1)*(V2-V1); %linear interpolation of basis

%interpolate full displacement vector at equilibrium
u_eq1 = ROMs.models{ind_p1}.thermal_eq;
u_eq2 = ROMs.models{ind_p2}.thermal_eq;
u_eq = u_eq1+(p-p1)/(p2-p1)*(u_eq2-u_eq1);

%% ROB
%use constant reduced order basis (corresponding to the VMs and modal
%derivatives of the cold structure
V =  ROMs.models{1}.V;

%% compute ROM
fullAssembly = ROMs.fullAssembly;
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

%norm of the residual to end Newton's iterations
c0 = norm(F_inertial) + norm(F_damping) + norm(F_int_r) + norm(F_ext_r);
end