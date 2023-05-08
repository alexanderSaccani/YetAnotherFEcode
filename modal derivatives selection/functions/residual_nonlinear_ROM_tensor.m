function [ r, drdqdd,drdqd,drdq, c0] = residual_nonlinear_ROM_tensor( q, qd, qdd, t, RedAssembly, tensors, Fext)

% Retrieve Mass and Damping matrices from reduced assembly
Mr = RedAssembly.DATA.M;
Cr = RedAssembly.DATA.C;

% Calculate nonlinear internal force vector and its derivatives
K2t = double(tensors.Q2); 
K3qq = double(ttv(tensors.Q3, {q, q}, [3,2])); 
K4qqq = double(ttv(tensors.Q4, {q, q, q}, [4,3,2]));
K3tq = double(ttv(tensors.Q3t, q, 3));
K4tqq = double(ttv(tensors.Q4t, {q,q}, [4,3]));

F_int = K2t * q + K3qq + K4qqq;
F_inertial = Mr * qdd;
F_damping = Cr * qd;
F_ext = Fext(t); %  External force vector at time t

% Calculate residual and its derivatives wrt qdd, qd, q
r = F_inertial +  F_damping  + F_int - F_ext ;
drdqdd = Mr;
drdqd =  Cr;
drdq = K2t + K3tq + K4tqq; 

% Calculate norm of the force components for convergence criteria
c0 = norm(F_inertial) + norm(F_damping) + norm(F_int) + norm(F_ext);

end