function [ r, drdqdd,drdqd,drdq, c0] = residual_thermal_reduced_nonlinear( q, qd, qdd, t, redAssembly, Fext, varargin)
%varargin should contain in order T, gradT (if required by the element)

V = redAssembly.V; %unconstrained basis 
M_V = redAssembly.DATA.M;
C_V = redAssembly.DATA.C;


if ~isempty(varargin)
    Tfun = varargin{1};
    varargTanStiff{1}  = Tfun(t);
end

if length(varargin) > 1
    gradTfun = varargin{2};
    varargTanStiff{2} = gradTfun(t);
end


u = V*q; %unconstrained vector of displacements
[K_V, F_V] = redAssembly.tangent_stiffness_and_force(u,varargTanStiff{:});

F_inertial = M_V * qdd;
F_damping = C_V * qd;
F_ext_V =  V.'*Fext(t);
r = F_inertial + F_damping + F_V - F_ext_V ;

drdqdd = M_V;
drdqd = C_V;
drdq = K_V;

c0 = norm(F_inertial) + norm(F_damping) + norm(F_V) + norm(F_ext_V);
end