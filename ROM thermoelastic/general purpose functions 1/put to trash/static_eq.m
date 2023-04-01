function [ u_lin, u ] = static_eq( Assembly, Fext, varargin )
% finds the equilibrium configuration of the model subject to Fext load and T.

[nsteps,tol,method,displayoption,vararginTanStiffForce] = parse_inputs(varargin{:});

K = Assembly.DATA.K;
F0 = Assembly.DATA.F0;

u_lin = Assembly.solve_system(K,Fext,F0); % in the linearized model, the internal forces are Kx + F0
u0 = Assembly.constrain_vector(u_lin);

switch method
    case 'fsolve'
        options = optimoptions('fsolve','SpecifyObjectiveGradient',true,...
            'MaxIterations',10000,'Display',displayoption);
        [ueq] = fsolve(@(u)f(u,Assembly,Fext,vararginTanStiffForce),u0,options);
        u = Assembly.unconstrain_vector(ueq);
        
    case 'newton'
        u = u_guess/nsteps;
        figure; xlabel('Normalized load');ylabel('$$\|\mathbf{u}\|$$','Interpreter','latex')
        h = animatedline;
        addpoints(h,0,0);
        for j = 1:nsteps
            Fext_j = j*Fext/nsteps;
            c0 = norm(Assembly.constrain_vector(Fext_j));
            it = 0;
            while true
                if isempty(vararginTanStiffForce)
                    [K, Fint] = Assembly.tangent_stiffness_and_force(x);
                else
                    [K, Fint] = Assembly.tangent_stiffness_and_force(x,vararginTanStiffForce{:});
                end
                residual = Fext_j - Fint;
                c = norm(Assembly.constrain_vector(residual))/c0;
                fprintf('STEP %d, ITERATION %d, RESIDUAL %d \n',j,it,c);
                if c < tol
                    break
                end
                correction = Assembly.solve_system(K,residual);
                u = u + correction;
                it = it + 1;
            end
            addpoints(h,j/nsteps,norm(u));
            drawnow 
        end
end

end

function [F,K] = f(u,Assembly,Fext,vararginTanStiffForce)
x = Assembly.unconstrain_vector(u);
    if isempty(vararginTanStiffForce)
        [Kt, Fint] = Assembly.tangent_stiffness_and_force(x);
    else
        [Kt, Fint] = Assembly.tangent_stiffness_and_force(x,vararginTanStiffForce{:});
    end
K = Assembly.constrain_matrix(Kt);
F = Assembly.constrain_vector(Fint - Fext);
end


function [nsteps,tol,method,displayoption,vararginTanStiffForce] = parse_inputs(varargin)
%% parsing inputs
defaultnsteps = 100;
defaulttol = 1e-6;
defaultmethod = 'fsolve';
defaultdisplay = 'final';
defaultVararginTanStiffForce = [];

p = inputParser;
addParameter(p,'nsteps',defaultnsteps, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'tol',defaulttol, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','positive'}) );
addParameter(p,'method',defaultmethod,@(x)validateattributes(x, ...
    {'char'},{'nonempty'}));
addParameter(p,'display',defaultdisplay,@(x)validateattributes(x, ...
    {'char'},{'nonempty'}));
addParameter(p,'vararginTanStiffForce', defaultVararginTanStiffForce);
parse(p,varargin{:});

nsteps = p.Results.nsteps;
tol = p.Results.tol;
method = p.Results.method;
displayoption = p.Results.display;
vararginTanStiffForce = p.Results.vararginTanStiffForce;

end