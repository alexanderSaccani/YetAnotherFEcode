function [ u_lin, u ] = static_equilibrium_thermal( Assembly, Fext, T, varargin )
% finds the equilibrium configuration of the model subject to Fext load.

[nsteps,tol,method,displayoption,initialGuess] = parse_inputs(varargin{:});

K = Assembly.DATA.K;
u_lin = Assembly.solve_system(K,Fext);

if isempty(initialGuess)
    u_guess = u_lin;
    u0 = Assembly.constrain_vector(u_lin);
else
    u_guess = initialGuess;  %unconstrained vector (called lin, but not necesserely it is the linear solution)
    u0 = Assembly.constrain_vector(u_guess); %constrained vector
end

switch method
    case 'fsolve'
        options = optimoptions('fsolve','SpecifyObjectiveGradient',true,...
            'MaxIterations',10000,'Display',displayoption);
        [ueq] = fsolve(@(u)f(u,Assembly,Fext,T),u0,options);
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
                [K, Fint] = Assembly.tangent_stiffness_and_force(u,T);
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

function [F,K] = f(u,Assembly,Fext,T)
x = Assembly.unconstrain_vector(u);
[Kt, Fint] = Assembly.tangent_stiffness_and_force(x,T);
K = Assembly.constrain_matrix(Kt);
F = Assembly.constrain_vector(Fint - Fext);
end


function [nsteps,tol,method,displayoption,initialGuess] = parse_inputs(varargin)
%% parsing inputs
defaultnsteps = 100;
defaulttol = 1e-6;
defaultmethod = 'fsolve';
defaultdisplay = 'final';
defaultInitialGuess = [];

p = inputParser;
addParameter(p,'nsteps',defaultnsteps, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'tol',defaulttol, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','positive'}) );
addParameter(p,'method',defaultmethod,@(x)validateattributes(x, ...
    {'char'},{'nonempty'}))
addParameter(p,'display',defaultdisplay,@(x)validateattributes(x, ...
    {'char'},{'nonempty'}))
addParameter(p,'initialGuess', defaultInitialGuess, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty'}) );
parse(p,varargin{:});

nsteps = p.Results.nsteps;
tol = p.Results.tol;
method = p.Results.method;
displayoption = p.Results.display;
initialGuess = p.Results.initialGuess;
end