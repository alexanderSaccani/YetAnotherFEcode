classdef ImplicitNewmarkRed < handle
    % Class for Implicit Newmark time Integration with reduced model
    % cf. flowchart 7.21 in the following reference:
    % Geradin, M., & Rixen, D. (2015). Mechanical Vibrations : Theory and Application to Structural Dynamics (Third Edition). 
    % Variable basis approach
    % Author: Alexander Saccani, Phd Candidate ETHz, May 2023
    properties
        % Time integration parameters
        alpha = 0.005
        beta
        h = 0 % time step size
        gamma
        tol = 1e-6      % Relative error tolerance
        
        
        Solution        % Solution data structure
        MaxNRit = 10
        ATS = false     % where adaptive time stepping should be on (true) or not (false)
        hmin = 0        % minimum timestep size (only used when ATS = true)
        NROpt = 3       % Maximum no. of N-R Iterations
        linear = false  % whether system is linear or not
    end
    
    methods
        function TI = ImplicitNewmarkRed(varargin)
            %Parsint the input
            p = inputParser;
            addParameter(p,'timestep',TI.h, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
            addParameter(p,'alpha',TI.alpha, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
            addParameter(p,'RelTol',TI.tol, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
            addParameter(p,'MaxNRit',TI.MaxNRit, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );
            addParameter(p,'linear', TI.linear, @(x)validateattributes(x,{'logical'},{'nonempty'}));
            addParameter(p,'hmin', TI.hmin, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
            addParameter(p,'ATS', TI.ATS, @(x)validateattributes(x,{'logical'},{'nonempty'}));
            
            parse(p,varargin{:});
            
            %Property assignements
            TI.alpha = p.Results.alpha;
            TI.beta = (1 + TI.alpha)^2/4;
            TI.gamma = 1/2 + TI.alpha;
            
            TI.h = p.Results.timestep;
            TI.tol = p.Results.RelTol;
            TI.MaxNRit = p.Results.MaxNRit;
            TI.hmin = p.Results.hmin;
            TI.linear = p.Results.linear;
        end
        
        function Integrate(obj,x0,xd0,xdd0,tmax, Residual,Vfun)            
            % Integrates with Initial condition x0,xd0 from [0 tmax]
            % Residual is a function handle that has the following syntax

            if obj.h ==0
                error('Please specify a positive time step')
            end
            
            tic
            t=0;
            time = t;
            
            %initial conditions
            x_old = x0;
            xd_old = xd0;
            xdd_old = xdd0;
            
            %save initial conditions in solution
            x = x0;
            xd = xd0;
            xdd = xdd0;
            
            NR = 0;
            R = 0;
            i = 1;    
            
            %initialize basis ID
            [~,basisIDold] = Vfun(t);
            
            while t < tmax
                t = t+obj.h;
                i = i+1;

                %load basis
                [V,basisIDnew] = Vfun(t);
                dimROM = size(V,2);
                
                %project the initial velocity and acceleration on the
                %current basis, if the basis changes. 
                if basisIDold ~= basisIDnew
                    xd_old = V*V'*xd_old;
                    xdd_old = V*V'*xdd_old;
                end       
                
                %initial condition in delta displacement, vel and acc are 0
%                 del_q_new = zeros(dimROM,1);
%                 del_qd_new = del_q_new;
                del_qdd_new =  zeros(dimROM,1);
                           
                it = -1; % iteration counter
                 
                %Newton Raphson iterations (in this loop we start from the
                %converged solution of the previous time step
                %(xold,xdold,xddold) and update it through newton rapshon
                %iterations to obtain the new solution at the current time
                %step (xnew,xdnew,xddnew)
                while true 
                    
                    it = it + 1;
                    
                    %update full order displacements, velocities and
                    %accelerations with the new delta acceleration in
                    %reduced space del_qdd_new
                    x_new = x_old + obj.h*xd_old + 0.5*obj.h^2*(xdd_old+2*obj.beta*V*del_qdd_new); 
                    xd_new = xd_old + obj.h*(xdd_old + obj.gamma*V*del_qdd_new);
                    xdd_new = xdd_old + V*del_qdd_new;
                    
                    %Evaluate residual (solve for del_qdd)
                    [r, drduddV, drdudV, drduV, c0] = Residual(x_new,xd_new,xdd_new,t,V);                        

                    %Check convergence
                    epsilon = norm(r)/c0;
                    disp(['Iteration ' num2str(it) ', Residual norm = '  num2str(epsilon)])
            
                    if (epsilon<obj.tol)  %Error < Tolerance : break
                        break;
                    else % Error >= Tolerance : perform correction
                        S = drduddV + obj.gamma * obj.h * drdudV + obj.beta * obj.h^2 * drduV;
                        D_del_qdd = -S\r; %correction of del_qdd acceleration 
                        del_qdd_new = del_qdd_new + D_del_qdd;
                    end

                    %Adapt time step to maintain an optimal number (obj.NROpt) of N-R iterations 
                    if obj.h > obj.hmin && obj.ATS
                        obj.h = max(obj.h*obj.NROpt/it, obj.hmin);
                    end

                    % When too many iterations
                    if (it > obj.MaxNRit)
                        warning('Max N-R iterations reached')                            
                        if  epsilon > 1 
                            disp('Exiting time integration: Too high a residual')
                            soltime=toc;
                            obj.Solution.time = time;
                            obj.Solution.q = x;
                            obj.Solution.qd = xd;
                            obj.Solution.qdd = xdd;
                            obj.Solution.NR = NR;
                            obj.Solution.R = R;
                            obj.Solution.soltime = soltime;
                            return
                        else
                            disp('Continuing iterations anyway since the residual is low')                            
                        end
                    end                    

                end
                
                % Update solution
                time = [time t];
                NR = [NR it];
                R = [R epsilon];
                disp(['time integration completed: ', num2str(100* t/tmax), '%'])
                
                x = [x x_new];
                xd = [xd xd_new];
                xdd = [xdd xdd_new]; 
                
                x_old = x_new;
                xd_old = xd_new;
                xdd_old = xdd_new;
                
            end
            
            soltime = toc;
            obj.Solution.time = time;
            obj.Solution.q = x;
            obj.Solution.qd = xd;
            obj.Solution.qdd = xdd;
            obj.Solution.NR = NR;
            obj.Solution.R = R;
            obj.Solution.soltime = soltime;
            
        end
        
%         function[q,qd,qdd] = Prediction(obj,q0,qd0,qdd0)
%             qd = qd0 + obj.h * (1 - obj.gamma) * qdd0;
%             q = q0 + obj.h * qd0 + (0.5-obj.beta) * obj.h^2 * qdd0;
%             qdd = zeros(length(q0),1);
%         end

%         function [q,qd,qdd] = Correction(obj,q,qd,qdd,Da)
%             q = q + obj.beta * obj.h^2 * Da;
%             qd = qd + obj.gamma * obj.h * Da;
%             qdd = qdd + Da; 
%         end
    end
end
