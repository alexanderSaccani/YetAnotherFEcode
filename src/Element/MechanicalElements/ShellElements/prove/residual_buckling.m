
function [r,drdX] = residual_buckling(X,Assembly,F0,varargin)
    % varargin -> put the extra arguments required (in addition to
    % displacement u) to the computation of the tangent stiffness and force
    % at element level

    u = X(1:end-1); %physical displacements (constrained)
    lam = X(end); %force multiplier (continuation parameter)
    
    u_f = Assembly.unconstrain_vector(u);
    [Kt, Fint] = Assembly.tangent_stiffness_and_force(u_f,varargin{:});
    
    %compute the residual
    r = Assembly.constrain_vector(Fint - lam*F0);
    
    %compute the derivative of the residual w.r.t. X
    K_c = Assembly.constrain_matrix(Kt);
    F0_c = Assembly.constrain_vector(F0);
    
    %output the derivative of the residual
    drdX = [K_c,-F0_c];
    
end