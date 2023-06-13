function [ r, drdqdd,drdqd,drdq, c0] = thermal_residual_interpolation( q, qd, qdd, t, Assembly, ROMs, pfunt, Fext, varargin)

%varargin{1} -> T(t)
%varargin{2} -> gradT(t)

Tfun = varargin{1};
varArgTanStiff{1} = Tfun(t);

if length(varargin) > 1
    gradTfun = varargin{2};
    varArgTanStiff{2} = gradTfun(t);
end
    
%interpolate the basis
p_current = pfunt(t);
xc = p_current(1);
psamples = ROMs.psamples;

%compare 
%p is column vector
nSamples = size(psamples,2);
dist_p_psamples1 = sqrt(sum((psamples - repmat(p_current,1,nSamples)).^2));
[dist1,closestModel1] = min(dist_p_psamples1);
V1 = ROMs.models{closestModel1}.V;
p1 = psamples(:,closestModel1);
x1 = p1(1);

disp(['model used is: ',num2str(closestModel1)])

psamples(:,closestModel1) = psamples(:,closestModel1)*inf;
dist_p_psamples2 = sqrt(sum((psamples - repmat(p_current,1,nSamples)).^2));
[dist2,closestModel2] = min(dist_p_psamples2);
V2 = ROMs.models{closestModel2}.V;
p2 = psamples(:,closestModel2);
x2 = p2(1);

disp(['model used is: ',num2str(closestModel2)])

distp1p2 = sqrt(sum((p1-p2).^2));

%interpolate the basis
%V = dist2/distp1p2*V1 + dist1/distp1p2*V2;
V = V1 + (xc-x1)/(x2-x1)*(V2-V1);



%construct the residual with the interpolated basis
%constrain the basis

M = Assembly.DATA.M;
C = Assembly.DATA.C;

u = V*q;
[K, F] = Assembly.tangent_stiffness_and_force(u,varArgTanStiff{:}); 


K_red = V'*K*V;
M_red = V'*M*V;
C_red = V'*C*V;
F_elastic = V'*F;
F_external =  V'*Fext(t); 

F_inertial = M_red * qdd;
F_damping = C_red * qd;
r = F_inertial + F_damping + F_elastic - F_external ;
drdqdd = M_red;
drdqd = C_red;
drdq = K_red;

c0 = norm(F_inertial) + norm(F_damping) + norm(F_elastic) + norm(F_external);

end