function [ r, drdqdd, drdqd, drdq] = residual_thermal_linear(q,qd,qdd,t,Assembly,Fext,varargin)

M = Assembly.DATA.M;
C = Assembly.DATA.C;

%the stiffness of the linearized model, is the tangent stiffness at
%equilibrium (null displacements, T varying). The thermal internal load vector
% is the thermal internal force for null displacements.
Tfun = varargin{1};
varArgTanStiffForc{1} = Tfun(t);
if length(varargin)>1
    gradTfun = varargin{2};
    varArgTanStiffForc{2} = gradTfun(t);
end

nDofsF = size(M,1);
[K_lin, gth] = Assembly.tangent_stiffness_and_force(zeros(nDofsF,1),varArgTanStiffForc{:});

% These matrices and the external forcing vector are appropriately constrained 
% according to the boundary conditions:
M_red = Assembly.constrain_matrix(M);
C_red = Assembly.constrain_matrix(C);
K_red = Assembly.constrain_matrix(K_lin);
F_red = Assembly.constrain_vector(Fext(t)-gth);

% Residual is computed according to the formula above:
r = M_red * qdd + C_red * qd + K_red * q - F_red ;
drdqdd = M_red;
drdqd = C_red;
drdq = K_red;

end