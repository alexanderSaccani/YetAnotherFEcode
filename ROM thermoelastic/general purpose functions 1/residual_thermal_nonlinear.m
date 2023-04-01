function [ r, drdqdd,drdqd,drdq, c0] = residual_thermal_nonlinear( q, qd, qdd, t, Assembly, Fext, varargin)

%varargin{1} -> T(t)
%varargin{2} -> gradT(t)

Tfun = varargin{1};
varArgTanStiff{1} = Tfun(t);

if length(varargin) > 1
    gradTfun = varargin{2};
    varArgTanStiff{2} = gradTfun(t);
end
    

M = Assembly.DATA.M;
C = Assembly.DATA.C;

u = Assembly.unconstrain_vector(q);
[K, F] = Assembly.tangent_stiffness_and_force(u,varArgTanStiff{:}); 

M_red = Assembly.constrain_matrix(M);
C_red = Assembly.constrain_matrix(C);
K_red = Assembly.constrain_matrix(K);
F_elastic = Assembly.constrain_vector(F);
F_external =  Assembly.constrain_vector(Fext(t)); 

F_inertial = M_red * qdd;
F_damping = C_red * qd;
r = F_inertial + F_damping + F_elastic - F_external ;
drdqdd = M_red;
drdqd = C_red;
drdq = K_red;

c0 = norm(F_inertial) + norm(F_damping) + norm(F_elastic) + norm(F_external);

end