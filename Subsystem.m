classdef Subsystem < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        index
        revisedIndex
        location
        neighbors
        exclusiveNeighbors
        outNeighbors 
        exclusiveOutNeighbors
        leastNeighbor  %\Gamma_i
        leastOutNeighbor %\Delta_i
        
        graphicHandles
        
        % stateSpaceDimentions
        dim_n % of state x
        dim_p % of input u
        dim_q % of disturbance w
        dim_m % of output y
        
        % Subsystem state-space variables
        x
        u
        w
        y
        newStateVariables
        
        % Subsystem state-space parameters
        % (these are subsystem-specific parameters, like : A = [A_{ij}:j\in\N_N]
        A = {}
        B = {}
        C = {}
        D = {}
        E = {}
        F = {}
        
        % controller gains
        localSFBLQRControllerGains % local state feedback controller gains
    end
    
    methods
        
        function obj = Subsystem(index,location)
            obj.index = index;
            obj.revisedIndex = index; % for now
            obj.location = location;
            
            % stateSpaceDimentions
            obj.dim_n = 4; % x
            obj.dim_p = 2; % u
            obj.dim_q = 2; % w
            obj.dim_m = 2; % y
        end
        
        function obj = loadParameters(obj,subsystems)
            n = obj.dim_n; % x
            p = obj.dim_p; % u
            q = obj.dim_q; % w
            m = obj.dim_m; % y
            % rss(n,p,m,2,numOfSubsystems) % random stable state spaces
            
            for j = 1:1:length(subsystems)
                if sum(j == obj.neighbors)==1 % j is a neighbor!
                    if j==obj.index
                        obj.A{j} = 10*(rand(n,n)-0.5); % A_ii
                        obj.B{j} = 7*(rand(n,p)-0.5);
                        obj.C{j} = 5*(rand(m,n)-0.5);
                        obj.D{j} = 2*(rand(m,p)-0.5);
                        obj.E{j} = 1*(rand(n,q)-0.5);
                        obj.F{j} = 0.1*(rand(m,q)-0.5);
                    else
                        [n_j, p_j, q_j, m_j] = subsystems(j).getDimensions();
                        obj.A{j} = 5*(rand(n,n_j)-0.5); % A_ij
                        obj.B{j} = 3*(rand(n,p_j)-0.5); 
                        obj.C{j} = 2*(rand(m,n_j)-0.5);
                        obj.D{j} = 1*(rand(m,p_j)-0.5);
                        obj.E{j} = 0.5*(rand(n,q_j)-0.5);
                        obj.F{j} = 0.05*(rand(m,q_j)-0.5);
                    end
                else
                    [n_j, p_j, q_j, m_j] = subsystems(j).getDimensions();
                    obj.A{j} = zeros(n,n_j); % A_ij
                    obj.B{j} = zeros(n,p_j); 
                    obj.C{j} = zeros(m,n_j);
                    obj.D{j} = zeros(m,p_j);
                    obj.E{j} = zeros(n,q_j);
                    obj.F{j} = zeros(m,q_j);
                end
            end
            
            obj.x = 10*rand(n,1); % initial state            
            obj.u = zeros(obj.dim_p,1); % initial control
            obj.w = zeros(obj.dim_q,1); % initial disturbance
        end
        
        
        function obj = loadStableParameters(obj,subsystems)
            n = obj.dim_n; % x
            p = obj.dim_p; % u
            q = obj.dim_q; % w
            m = obj.dim_m; % y
            
%             rss(n,p,m,2,length(subsystems)) % random stable state spaces
            
            for j = 1:1:length(subsystems)
                if sum(j == obj.neighbors)==1 % j is a neighbor!
                    if j==obj.index
                        sys1 = rss(n,p,m);
                        sys2 = rss(n,q,m);
                        obj.A{j} = sys1.A; % A_ii
                        obj.B{j} = sys1.B;
                        obj.C{j} = sys1.C;
                        obj.D{j} = sys1.D;
                        obj.E{j} = 0.5*sys2.B;
                        obj.F{j} = 0.1*sys2.D;
                    else
                        [n_j, p_j, q_j, m_j] = subsystems(j).getDimensions();
                        sys1 = rss(n,p_j,m_j);  % assume n=n_j
                        sys2 = rss(n,q_j,m_j) ; % assume n=n_j
                        obj.A{j} = 0.1*sys1.A; % A_ij
                        obj.B{j} = 0.1*sys1.B; 
                        obj.C{j} = 0.1*sys1.C;
                        obj.D{j} = 0.1*sys1.D;
                        obj.E{j} = 0.05*sys2.B;
                        obj.F{j} = 0.01*sys2.D;
                    end
                else
                    [n_j, p_j, q_j, m_j] = subsystems(j).getDimensions();
                    obj.A{j} = zeros(n,n_j); % A_ij
                    obj.B{j} = zeros(n,p_j); 
                    obj.C{j} = zeros(m,n_j);
                    obj.D{j} = zeros(m,p_j);
                    obj.E{j} = zeros(n,q_j);
                    obj.F{j} = zeros(m,q_j);
                end
            end
            
            obj.x = 1*rand(n,1); % initial state            
            obj.u = zeros(obj.dim_p,1); % initial control
            obj.w = zeros(obj.dim_q,1); % initial disturbance
        end
        
        function output = designLocalSFBLQRControllerGains(obj)
            A = obj.A{obj.index};
            B = obj.B{obj.index};
            Q = 100*eye(obj.dim_n);
            R = 20*eye(obj.dim_p);
            
            [K,~,~] = lqr(A,B,Q,R);
            obj.localSFBLQRControllerGains = K;
        end
        
        function [n,p,q,m] = getDimensions(obj) %size of state, input, disturbance and output
            n = obj.dim_n;
            p = obj.dim_p;
            q = obj.dim_q;
            m = obj.dim_m;
        end
        
        function [x,u,w,y] = getStateVariables(obj) % state x_i , input u_i , disturbance w_i , output y_i
            x = obj.x;
            u = obj.u;
            w = obj.w;
            y = obj.y;
        end
        
        function outputArg = drawSubsystem(obj)
            hold on
            viscircles(obj.location,0.015,'Color','b');
            text(obj.location(1)-0.02,obj.location(2)-0.04,num2str(obj.index),'Color','k','FontSize',10);
        end
        
        function outputArg = drawIndex(obj,indexVal)
            hold on
            text(obj.location(1)-0.04,obj.location(2),num2str(indexVal),'Color','r','FontSize',10);
        end 
        
        
        function [xCost, uCost, data] = update(obj,deltaT,subsystems)
            
            % find new x based on current x, u and w
                         
            % Computing sumAx sumBx sumEx sumCx sumDx sumFx
            sumAx = zeros(obj.dim_n,1);
            sumBu = zeros(obj.dim_n,1);
            sumCx = zeros(obj.dim_m,1);
            sumDu = zeros(obj.dim_m,1);
            sumEw = zeros(obj.dim_n,1);
            sumFw = zeros(obj.dim_m,1);
            
            for jInd = 1:1:length(obj.neighbors)
                j = obj.neighbors(jInd); % j=i is also a possibility
                sumAx = sumAx + obj.A{j}*subsystems(j).x;
                sumBu = sumBu + obj.B{j}*subsystems(j).u;
                sumCx = sumCx + obj.C{j}*subsystems(j).x;
                sumDu = sumDu + obj.D{j}*subsystems(j).u;
                sumEw = sumEw + obj.E{j}*subsystems(j).w;
                sumFw = sumFw + obj.F{j}*subsystems(j).w;
            end
            
            % next state and output
            xNew = obj.x + deltaT*(sumAx+sumBu+sumEw);
            yNew = sumCx+sumDu+sumFw;
           
            % control u and disturbance w for the next iteration
            uNew = -obj.localSFBLQRControllerGains*xNew; % zeros(obj.dim_p,1);
            wNew = 0.05*randn(obj.dim_q,1);
            
            % compiling new data
            xCost = deltaT*xNew'*xNew;
            uCost = deltaT*uNew'*uNew;
            data = [norm(xNew), norm(uNew), norm(wNew), norm(yNew)];
            
            
            obj.newStateVariables.xNew = xNew;
            obj.newStateVariables.yNew = yNew;
            obj.newStateVariables.uNew = uNew;
            obj.newStateVariables.wNew = wNew;
            
        end
        
        function outputArg = finishUpdate(obj)
            
            obj.x = obj.newStateVariables.xNew;
            obj.y = obj.newStateVariables.yNew;
            obj.u = obj.newStateVariables.uNew;
            obj.w = obj.newStateVariables.wNew;
            
        end
    end
end

