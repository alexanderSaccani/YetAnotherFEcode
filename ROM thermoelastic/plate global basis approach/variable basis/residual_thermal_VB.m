function [ r, drdqdd,drdqd,drdq, c0] = residual_thermal_VB( x, xd, xdd, t, Assembly, V, W, Fext, varargin)
%varargin should contain in order T, gradT (if required by the element)

%q, qd, qdd are full displacements unconstrained

%unconstrained basis 
M = Assembly.DATA.M; 
C = Assembly.DATA.C;


if ~isempty(varargin)
    Tfun = varargin{1};
    varargTanStiff{1}  = Tfun(t);
end

if length(varargin) > 1
    gradTfun = varargin{2};
    varargTanStiff{2} = gradTfun(t);
end

[K, F_internal] = Assembly.tangent_stiffness_and_force(x,varargTanStiff{:});

F_inertial = M * xdd;
F_damping = C* xd;
F_ext =  Fext(t);

if W.flag == true
    L = W.basis;
else
    L = V;
end

r = L'*(F_inertial + F_damping + F_internal - F_ext) ;

drdqdd = L'*M*V; %can be put in input
drdqd = L'*C*V; %can be put in input
drdq = L'*K*V;

c0 = norm(L'*F_inertial) + norm(L'*F_internal) + norm(L'*F_ext) +  norm(L'*F_damping);

end