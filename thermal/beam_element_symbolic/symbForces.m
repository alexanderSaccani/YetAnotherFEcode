% Alexander Saccani - PhD Candidate ETH Zurich - 2022
% symbolic generation of shape functions and internal forces for thermal
% beam element
clearvars
clc

syms xi;
syms u1 v1 w1 u2 v2 w2 T1 T2;
syms l A E I alpha;

p = [u1,v1,w1,u2,v2,w2].'; %nodal displacements

%% define shape functions__________________________________________________

%shape functions for horizontal displacements
hu = [1/2-xi/2, 0, 0, 1/2+xi/2, 0, 0];

bu = diff(hu,xi);
% subs(hu,xi,-1)
% subs(hu,xi,+1)

%shape functions for vertical displacements
B = [1 -1 1 -1; 1 1 1 1; 0 2/l -4/l 6/l; 0 2/l 4/l 6/l];
B_inv = inv(B);

hv = zeros(1,6);

for ii = 0:3

    dhv = [0, B_inv(ii+1,1)*xi^ii, B_inv(ii+1,3)*xi^ii, 0, B_inv(ii+1,2)*xi^ii, B_inv(ii+1,4)*xi^ii ];
    hv = hv + dhv;

end

%first order derivative
bv = diff(hv,xi);

%second order derivative
cv = diff(bv,xi);

% subs(hv,xi,-1)
% subs(hv,xi,+1)
% subs(bv,xi,-1)
% subs(bv,xi,+1)

%shape function for temperature
hT = [1/2-xi/2, 1/2+xi/2];

%% define tensors of internal forces_______________________________________

%tensors of terms linear in displacements
L = zeros(6,6);
L = sym(L);

xi_int = [-1,1]; %integration interval

for i = 1:6
    for j = 1:6
    L(i,j) = int( 2*A*E/l*bu(i)*bu(j), xi, xi_int) + ...
        int( 8*E*I/l^3*cv(i)*cv(j), xi, xi_int);
    end
end

LT = zeros(6,6,2);
LT = sym(LT);
for i = 1:6
    for j = 1:6
        for k = 1:2
           LT(i,j,k) = int( -2*A*E*alpha/l*bv(i)*bv(j)*hT(k), xi, xi_int);
        end
    end
end

Q = zeros(6,6,6);
Q = sym(Q);
for i = 1:6
    for j = 1:6
        for k = 1:6
            Q(i,j,k) = int( 4*A*E/l^2*bv(i)*bv(k)*bu(j), xi, xi_int) + ...
                int( 2*A*E/l^2*bu(i)*bv(j)*bv(k), xi, xi_int);
        end
    end
end

C = zeros(6,6,6,6);
C = sym(C);
for i = 1:6
    for j = 1:6
        for k = 1:6
            for m = 1:6
               a = bv(i)*bv(j)*bv(k)*bv(m);
               C(i,j,k,m) = int( a*(4*A*E/l^3), xi, xi_int);
            end
        end
    end
end

gT = zeros(6,2);
gT = sym(gT);

for i = 1:6
    for j = 1:2
        gT(i,j) = int(-A*E*alpha*bu(i)*hT(j), xi, xi_int);
    end
end

%% assemble force vector

F = zeros(6,1);
F = sym(F);

F1 = F;
F2 = F;
F3 = F;
g = F;

for i = 1:6
    for j = 1:6
        F1(i) = F1(i) + L(i,j)*p(j) + LT(i,j,1)*T1*p(j) + LT(i,j,2)*T2*p(j);
    end
end

for i = 1:6
    for j = 1:6
        for k = 1:6
            F2(i) = F2(i) + Q(i,j,k)*p(j)*p(k);
        end
    end
end

for i = 1:6
    for j = 1:6
        for k = 1:6
            for m = 1:6
                F3(i) = F3(i) + C(i,j,k,m)*p(j)*p(k)*p(m);
            end
        end
    end
end

g = gT(:,1)*T1 + gT(:,2)*T2;

F = F1 + F2 + F3 + g ; %total internal force

%uniform body force
hf = zeros(6,3);
hf = sym(hf);

hf(:,1) = int(hu*l/2,xi,xi_int).';
hf(:,2) = int(hv*l/2,xi,xi_int).';
hf(:,3) = int(bv,xi,xi_int).';


%% compare results obtained with this code, and the one from Shobhit's code
syms t1 t2 
F_old =  [(A*E*(- 2*t1^2 + t1*t2 - 2*t2^2 + 15*T1*alpha + 15*T2*alpha))/30 - ((A*E*(18*w1^2 - 36*w1*w2 + 18*w2^2))/30 - (A*E*l*(30*u1 - 30*u2 - 3*t1*w1 + 3*t1*w2 - 3*t2*w1 + 3*t2*w2))/30)/l^2;
                (E*(3360*I*w1 - 3360*I*w2 + 288*A*w1^3 - 288*A*w2^3 - A*l^3*t1^3 - A*l^3*t2^3 + 1680*I*l*t1 + 1680*I*l*t2 + 864*A*w1*w2^2 - 864*A*w1^2*w2 - 28*A*l^2*t1*u1 + 28*A*l^2*t1*u2 - 28*A*l^2*t2*u1 + 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 + 3*A*l^3*t1*t2^2 + 3*A*l^3*t1^2*t2 + 36*A*l^2*t1^2*w1 - 36*A*l^2*t1^2*w2 + 36*A*l^2*t2^2*w1 - 36*A*l^2*t2^2*w2 - 336*A*l*u1*w1 + 336*A*l*u1*w2 + 336*A*l*u2*w1 - 336*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 216*A*l*t2*w1*w2 - 28*A*T1*alpha*l^3*t2 - 28*A*T2*alpha*l^3*t1 - 168*A*T1*alpha*l^2*w1 + 168*A*T1*alpha*l^2*w2 - 168*A*T2*alpha*l^2*w1 + 168*A*T2*alpha*l^2*w2))/(280*l^3);
                (E*(5040*I*w1 - 5040*I*w2 + 108*A*w1^3 - 108*A*w2^3 + 24*A*l^3*t1^3 - 3*A*l^3*t2^3 + 3360*I*l*t1 + 1680*I*l*t2 + 324*A*w1*w2^2 - 324*A*w1^2*w2 - 112*A*l^2*t1*u1 + 112*A*l^2*t1*u2 + 28*A*l^2*t2*u1 - 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 6*A*l^3*t1*t2^2 - 9*A*l^3*t1^2*t2 - 9*A*l^2*t1^2*w1 + 9*A*l^2*t1^2*w2 + 9*A*l^2*t2^2*w1 - 9*A*l^2*t2^2*w2 - 84*A*l*u1*w1 + 84*A*l*u1*w2 + 84*A*l*u2*w1 - 84*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 84*A*T1*alpha*l^3*t1 + 14*A*T1*alpha*l^3*t2 - 28*A*T2*alpha*l^3*t1 + 14*A*T2*alpha*l^3*t2 - 84*A*T2*alpha*l^2*w1 + 84*A*T2*alpha*l^2*w2 + 18*A*l^2*t1*t2*w1 - 18*A*l^2*t1*t2*w2))/(840*l^2);
                ((A*E*(18*w1^2 - 36*w1*w2 + 18*w2^2))/30 - (A*E*l*(30*u1 - 30*u2 - 3*t1*w1 + 3*t1*w2 - 3*t2*w1 + 3*t2*w2))/30)/l^2 - (A*E*(- 2*t1^2 + t1*t2 - 2*t2^2 + 15*T1*alpha + 15*T2*alpha))/30;
                -(E*(3360*I*w1 - 3360*I*w2 + 288*A*w1^3 - 288*A*w2^3 - A*l^3*t1^3 - A*l^3*t2^3 + 1680*I*l*t1 + 1680*I*l*t2 + 864*A*w1*w2^2 - 864*A*w1^2*w2 - 28*A*l^2*t1*u1 + 28*A*l^2*t1*u2 - 28*A*l^2*t2*u1 + 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 + 3*A*l^3*t1*t2^2 + 3*A*l^3*t1^2*t2 + 36*A*l^2*t1^2*w1 - 36*A*l^2*t1^2*w2 + 36*A*l^2*t2^2*w1 - 36*A*l^2*t2^2*w2 - 336*A*l*u1*w1 + 336*A*l*u1*w2 + 336*A*l*u2*w1 - 336*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 216*A*l*t2*w1*w2 - 28*A*T1*alpha*l^3*t2 - 28*A*T2*alpha*l^3*t1 - 168*A*T1*alpha*l^2*w1 + 168*A*T1*alpha*l^2*w2 - 168*A*T2*alpha*l^2*w1 + 168*A*T2*alpha*l^2*w2))/(280*l^3);
                (E*(5040*I*w1 - 5040*I*w2 + 108*A*w1^3 - 108*A*w2^3 - 3*A*l^3*t1^3 + 24*A*l^3*t2^3 + 1680*I*l*t1 + 3360*I*l*t2 + 324*A*w1*w2^2 - 324*A*w1^2*w2 + 28*A*l^2*t1*u1 - 28*A*l^2*t1*u2 - 112*A*l^2*t2*u1 + 112*A*l^2*t2*u2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 - 9*A*l^3*t1*t2^2 + 6*A*l^3*t1^2*t2 + 9*A*l^2*t1^2*w1 - 9*A*l^2*t1^2*w2 - 9*A*l^2*t2^2*w1 + 9*A*l^2*t2^2*w2 - 84*A*l*u1*w1 + 84*A*l*u1*w2 + 84*A*l*u2*w1 - 84*A*l*u2*w2 - 216*A*l*t2*w1*w2 + 14*A*T1*alpha*l^3*t1 - 28*A*T1*alpha*l^3*t2 + 14*A*T2*alpha*l^3*t1 - 84*A*T2*alpha*l^3*t2 - 84*A*T1*alpha*l^2*w1 + 84*A*T1*alpha*l^2*w2 + 18*A*l^2*t1*t2*w1 - 18*A*l^2*t1*t2*w2))/(840*l^2)];

F_old = subs(F_old,w1,v1);
F_old = subs(F_old,w2,v2);
F_old = subs(F_old,t1,w1);
F_old = subs(F_old,t2,w2);


simplify(F-F_old)

%tangent stiffness matrix

Ltot = L + LT(:,:,1)*T1 + LT(:,:,2)*T2;

K_tan = Ltot;

for i = 1:6
    for j = 1:6
        for k = 1:6
            K_tan(i,j) = K_tan(i,j) + Q(i,j,k)*p(k);
            K_tan(i,j) = K_tan(i,j) + Q(i,k,j)*p(k);
            for m = 1:6
                K_tan(i,j) = K_tan(i,j) + C(i,j,k,m)*p(k)*p(m);
                K_tan(i,j) = K_tan(i,j) + C(i,k,j,m)*p(k)*p(m);
                K_tan(i,j) = K_tan(i,j) + C(i,k,m,j)*p(k)*p(m);
            end
        end
    end
end

        
K_old =  [                                         (A*E)/l,                                                                                                                                                                                    -((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                          - (A*E*(4*t1 - t2))/30 - (A*E*(3*w1 - 3*w2))/(30*l),                                        -(A*E)/l,                                                                                                                                                                                     ((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                            (A*E*(t1 - 4*t2))/30 - (A*E*(3*w1 - 3*w2))/(30*l);
                -(A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2),  (E*(3360*I + 864*A*w1^2 + 864*A*w2^2 + 36*A*l^2*t1^2 + 36*A*l^2*t2^2 - 336*A*l*u1 + 336*A*l*u2 - 1728*A*w1*w2 - 168*A*T1*alpha*l^2 - 168*A*T2*alpha*l^2 + 216*A*l*t1*w1 - 216*A*l*t1*w2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(280*l^3),                                               (E*(1680*I + 108*A*w1^2 + 108*A*w2^2 - 3*A*l^2*t1^2 + 3*A*l^2*t2^2 - 28*A*l*u1 + 28*A*l*u2 - 216*A*w1*w2 - 28*A*T2*alpha*l^2 + 6*A*l^2*t1*t2 + 72*A*l*t1*w1 - 72*A*l*t1*w2))/(280*l^2),  (A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2), -(E*(3360*I + 864*A*w1^2 + 864*A*w2^2 + 36*A*l^2*t1^2 + 36*A*l^2*t2^2 - 336*A*l*u1 + 336*A*l*u2 - 1728*A*w1*w2 - 168*A*T1*alpha*l^2 - 168*A*T2*alpha*l^2 + 216*A*l*t1*w1 - 216*A*l*t1*w2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(280*l^3),                                               (E*(1680*I + 108*A*w1^2 + 108*A*w2^2 + 3*A*l^2*t1^2 - 3*A*l^2*t2^2 - 28*A*l*u1 + 28*A*l*u2 - 216*A*w1*w2 - 28*A*T1*alpha*l^2 + 6*A*l^2*t1*t2 + 72*A*l*t2*w1 - 72*A*l*t2*w2))/(280*l^2);
                -(A*E*(3*w1 - 3*w2 + 4*l*t1 - l*t2))/(30*l),                                              (E*(5040*I + 324*A*w1^2 + 324*A*w2^2 - 9*A*l^2*t1^2 + 9*A*l^2*t2^2 - 84*A*l*u1 + 84*A*l*u2 - 648*A*w1*w2 - 84*A*T2*alpha*l^2 + 18*A*l^2*t1*t2 + 216*A*l*t1*w1 - 216*A*l*t1*w2))/(840*l^2), (E*(1680*I + 54*A*w1^2 + 54*A*w2^2 + 36*A*l^2*t1^2 + 3*A*l^2*t2^2 - 56*A*l*u1 + 56*A*l*u2 - 108*A*w1*w2 - 42*A*T1*alpha*l^2 - 14*A*T2*alpha*l^2 - 9*A*l^2*t1*t2 - 9*A*l*t1*w1 + 9*A*l*t1*w2 + 9*A*l*t2*w1 - 9*A*l*t2*w2))/(420*l),    (A*E*(3*w1 - 3*w2 + 4*l*t1 - l*t2))/(30*l),                                             -(E*(5040*I + 324*A*w1^2 + 324*A*w2^2 - 9*A*l^2*t1^2 + 9*A*l^2*t2^2 - 84*A*l*u1 + 84*A*l*u2 - 648*A*w1*w2 - 84*A*T2*alpha*l^2 + 18*A*l^2*t1*t2 + 216*A*l*t1*w1 - 216*A*l*t1*w2))/(840*l^2),                                   (E*(1680*I - 9*A*l^2*t1^2 - 9*A*l^2*t2^2 + 28*A*l*u1 - 28*A*l*u2 + 14*A*T1*alpha*l^2 + 14*A*T2*alpha*l^2 + 12*A*l^2*t1*t2 + 18*A*l*t1*w1 - 18*A*l*t1*w2 + 18*A*l*t2*w1 - 18*A*l*t2*w2))/(840*l);
                -(A*E)/l,                                                                                                                                                                                     ((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                            (A*E*(4*t1 - t2))/30 + (A*E*(3*w1 - 3*w2))/(30*l),                                         (A*E)/l,                                                                                                                                                                                    -((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                            (A*E*(3*w1 - 3*w2))/(30*l) - (A*E*(t1 - 4*t2))/30;
                (A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2), -(E*(3360*I + 864*A*w1^2 + 864*A*w2^2 + 36*A*l^2*t1^2 + 36*A*l^2*t2^2 - 336*A*l*u1 + 336*A*l*u2 - 1728*A*w1*w2 - 168*A*T1*alpha*l^2 - 168*A*T2*alpha*l^2 + 216*A*l*t1*w1 - 216*A*l*t1*w2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(280*l^3),                                              -(E*(1680*I + 108*A*w1^2 + 108*A*w2^2 - 3*A*l^2*t1^2 + 3*A*l^2*t2^2 - 28*A*l*u1 + 28*A*l*u2 - 216*A*w1*w2 - 28*A*T2*alpha*l^2 + 6*A*l^2*t1*t2 + 72*A*l*t1*w1 - 72*A*l*t1*w2))/(280*l^2), -(A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2),  (E*(3360*I + 864*A*w1^2 + 864*A*w2^2 + 36*A*l^2*t1^2 + 36*A*l^2*t2^2 - 336*A*l*u1 + 336*A*l*u2 - 1728*A*w1*w2 - 168*A*T1*alpha*l^2 - 168*A*T2*alpha*l^2 + 216*A*l*t1*w1 - 216*A*l*t1*w2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(280*l^3),                                              -(E*(1680*I + 108*A*w1^2 + 108*A*w2^2 + 3*A*l^2*t1^2 - 3*A*l^2*t2^2 - 28*A*l*u1 + 28*A*l*u2 - 216*A*w1*w2 - 28*A*T1*alpha*l^2 + 6*A*l^2*t1*t2 + 72*A*l*t2*w1 - 72*A*l*t2*w2))/(280*l^2);
                -(A*E*(3*w1 - 3*w2 - l*t1 + 4*l*t2))/(30*l),                                              (E*(5040*I + 324*A*w1^2 + 324*A*w2^2 + 9*A*l^2*t1^2 - 9*A*l^2*t2^2 - 84*A*l*u1 + 84*A*l*u2 - 648*A*w1*w2 - 84*A*T1*alpha*l^2 + 18*A*l^2*t1*t2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(840*l^2),                                   (E*(1680*I - 9*A*l^2*t1^2 - 9*A*l^2*t2^2 + 28*A*l*u1 - 28*A*l*u2 + 14*A*T1*alpha*l^2 + 14*A*T2*alpha*l^2 + 12*A*l^2*t1*t2 + 18*A*l*t1*w1 - 18*A*l*t1*w2 + 18*A*l*t2*w1 - 18*A*l*t2*w2))/(840*l),    (A*E*(3*w1 - 3*w2 - l*t1 + 4*l*t2))/(30*l),                                             -(E*(5040*I + 324*A*w1^2 + 324*A*w2^2 + 9*A*l^2*t1^2 - 9*A*l^2*t2^2 - 84*A*l*u1 + 84*A*l*u2 - 648*A*w1*w2 - 84*A*T1*alpha*l^2 + 18*A*l^2*t1*t2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(840*l^2), (E*(1680*I + 54*A*w1^2 + 54*A*w2^2 + 3*A*l^2*t1^2 + 36*A*l^2*t2^2 - 56*A*l*u1 + 56*A*l*u2 - 108*A*w1*w2 - 14*A*T1*alpha*l^2 - 42*A*T2*alpha*l^2 - 9*A*l^2*t1*t2 + 9*A*l*t1*w1 - 9*A*l*t1*w2 - 9*A*l*t2*w1 + 9*A*l*t2*w2))/(420*l)];


K_old = subs(K_old,w1,v1);
K_old = subs(K_old,w2,v2);
K_old = subs(K_old,t1,w1);
K_old = subs(K_old,t2,w2);


simplify(K_old - K_tan);






















