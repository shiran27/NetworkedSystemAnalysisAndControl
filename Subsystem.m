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
        dim_l % of performance metric z
        
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
        
        G = {}
        H = {}
        J = {}
        
        % controller gains
        controllerGains
        observerGains
        
        % data storage for distributed analysis and synthesis algorithms
        dataToBeDistributed % P_i and \tilde{W}_i in hte case of stability analysis
        
        % test network matrix
        testMatrix = {};
    end
    
    methods
        
        function obj = Subsystem(index,location,dims)
            obj.index = index;
            obj.revisedIndex = index; % for now
            obj.location = location;
            
            % stateSpaceDimentions
            obj.dim_n = dims.n; % x (n)
            obj.dim_p = dims.p; % u (p)
            obj.dim_q = dims.q; % w (q)
            obj.dim_m = dims.m; % y (m)
            obj.dim_l = dims.l; % z (l)
        end
        
        function obj = loadParameters(obj,subsystems)
            n = obj.dim_n; % x
            p = obj.dim_p; % u
            q = obj.dim_q; % w
            m = obj.dim_m; % y
            l = obj.dim_l; % z
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
                        [n_j, p_j, q_j, m_j, l_j] = subsystems(j).getDimensions();
                        obj.A{j} = 5*(rand(n,n_j)-0.5); % A_ij
                        obj.B{j} = 3*(rand(n,p_j)-0.5); 
                        obj.C{j} = 2*(rand(m,n_j)-0.5);
                        obj.D{j} = 1*(rand(m,p_j)-0.5);
                        obj.E{j} = 0.5*(rand(n,q_j)-0.5);
                        obj.F{j} = 0.05*(rand(m,q_j)-0.5);
                    end
                else
                    [n_j, p_j, q_j, m_j, l_j] = subsystems(j).getDimensions();
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
        
        
        function obj = loadStableParameters(obj,subsystems,diags)
            
            % diags = nature of [A,B,C,D,E,F,G,H,J], 1 diag, 0 general, -1 zero
            
            n = obj.dim_n; % x
            p = obj.dim_p; % u
            q = obj.dim_q; % w
            m = obj.dim_m; % y
            l = obj.dim_l; % z
            
%             rss(n,p,m,2,length(subsystems)) % random stable state spaces
            
            for j = 1:1:length(subsystems)
                if sum(j == obj.neighbors)==1 % j is a neighbor!
                    if j==obj.index % j = i
                        sys1 = rss(n,m,p); % n states, m outputs, p inputs 
                        sys2 = rss(n,m,q); % n states, m outputs, q inputs
                        obj.A{j} = sys1.A*(diags(1)>=0); % A_ii
                        obj.B{j} = sys1.B*(diags(2)>=0);
                        obj.C{j} = sys1.C*(diags(3)>=0);
                        obj.D{j} = sys1.D*(diags(4)>=0);
                        obj.E{j} = 0.01*sys2.B*(diags(5)>=0);
                        obj.F{j} = 0.01*sys2.D*(diags(6)>=0);
                        obj.G{j} = ones(l,n)*(diags(7)>=0);
                        obj.H{j} = ones(l,p)*(diags(8)>=0);
                        obj.J{j} = ones(l,q)*(diags(9)>=0);
                    else % j != i
                        [n_j, p_j, q_j, m_j, l_j] = subsystems(j).getDimensions();
                        sys1 = rss(n,m_j,p_j); % assume n=n_j states, m_j outputs, p_j inputs
                        sys2 = rss(n,m_j,q_j); % assume n=n_j states, m_j outputs, q_j inputs
                        obj.A{j} = 0.1*sys1.A*(diags(1)==0); % A_ij
                        obj.B{j} = 0.1*sys1.B*(diags(2)==0); 
                        obj.C{j} = 0.1*sys1.C*(diags(3)==0);
                        obj.D{j} = 0.1*sys1.D*(diags(4)==0);
                        obj.E{j} = 0.005*sys2.B*(diags(5)==0);
                        obj.F{j} = 0.005*sys2.D*(diags(6)==0);
                        obj.G{j} = -ones(l,n_j)*(diags(7)==0);
                        obj.H{j} = ones(l,p_j)*(diags(8)==0);
                        obj.J{j} = ones(l,q_j)*(diags(9)==0);
                    end
                else
                    [n_j, p_j, q_j, m_j, l_j] = subsystems(j).getDimensions();
                    obj.A{j} = zeros(n,n_j); % A_ij
                    obj.B{j} = zeros(n,p_j); 
                    obj.C{j} = zeros(m,n_j);
                    obj.D{j} = zeros(m,p_j);
                    obj.E{j} = zeros(n,q_j);
                    obj.F{j} = zeros(m,q_j);
                    obj.G{j} = zeros(l,n_j);
                    obj.H{j} = zeros(l,p_j);
                    obj.J{j} = zeros(l,q_j);
                end
            end
            
            obj.x = 10*rand(n,1); % initial state            
            obj.u = zeros(obj.dim_p,1); % initial control
            obj.w = zeros(obj.dim_q,1); % initial disturbance
        end
        
        
        
        function obj = loadPassiveParameters(obj,subsystems)
            n = obj.dim_n; % x
            p = obj.dim_p; % u
            q = obj.dim_q; % w
            m = obj.dim_m; % y
            
%             rss(n,p,m,2,length(subsystems)) % random stable state spaces
            
            for j = 1:1:length(subsystems)
                if sum(j == obj.neighbors)==1 % j is a neighbor!
                    if j==obj.index
                        sys1 = rss(n,m,p);
                        while ~isPassive(sys1)
                            sys1 = rss(n,m,p);
                        end
                        sys2 = rss(n,m,q);
                        obj.A{j} = sys1.A; % A_ii
                        obj.B{j} = 0.5*sys2.B;
                        obj.C{j} = sys1.C;
                        obj.D{j} = 0*sys2.D;
                        obj.E{j} = sys1.B;
                        obj.F{j} = sys1.D;
                    else
                        [n_j, p_j, q_j, m_j, l_j] = subsystems(j).getDimensions();
                        sys1 = rss(n,m_j,p_j);  % assume n=n_j
                        while ~isPassive(sys1)
                            sys1 = rss(n,m_j,p_j);  % assume n=n_j
                        end
                        sys2 = rss(n,m_j,q_j); % assume n=n_j
                        obj.A{j} = 0.1*sys1.A; % A_ij
                        obj.B{j} = 0*sys2.B; 
                        obj.C{j} = 0*sys1.C;
                        obj.D{j} = 0*sys2.D;
                        obj.E{j} = 0.1*sys1.B;
                        obj.F{j} = 0*sys1.D;
                    end
                else
                    [n_j, p_j, q_j, m_j, l_j] = subsystems(j).getDimensions();
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
            obj.controllerGains.localSFBLQR = K;
        end
        
        function [n,p,q,m,l] = getDimensions(obj) %size of state, input, disturbance and output
            n = obj.dim_n;
            p = obj.dim_p;
            q = obj.dim_q;
            m = obj.dim_m;
            l = obj.dim_l;
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
            text(obj.location(1)-0.035,obj.location(2)-0.03,num2str(obj.index),'Color','k','FontSize',10);
        end
        
        function outputArg = drawIndex(obj,indexVal,angleVal,colorVal)
            hold on
            r = 0.04;
            locX = obj.location(1)+r*cos(angleVal);
            locY = obj.location(2)+r*sin(angleVal);
            text(locX,locY,num2str(indexVal),'Color',colorVal,'FontSize',10);
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
            uNew = obj.getControlInput(xNew,subsystems,'uncontrolled')
            
            
            % disturbance w
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
        
        function uNew = getControlInput(obj,xNew,subsystems,method)
            if isequal(method,'uncontrolled')       % uncontrolled system
                uNew = zeros(obj.dim_p,1); 
            elseif isequal(method,'localSFBLQR')    % controlled based on local information A_ii, B_ii and SFB LQR
                uNew = -obj.controllerGains.localSFBLQR*xNew;
            elseif isequal(method,'gobalSFBLQR')
                uNew = zeros(obj.dim_p,1);          % control based on desgined global SFB LQR controller gains
                for jInd = 1:1:length(obj.neighbors)
                    j = obj.neighbors(jInd);
                    uNew = uNew - obj.controllerGains.globalSFBLQR{j}*subsystems(j).x;
                end
            elseif isequal(method,'gobalSFBStabLMI')
                
            elseif isequal(method,'gobalSFBDissLMI')
                
            elseif isequal(method,'localSFBStabLMI')
                
            elseif isequal(method,'localSFBDissLMI')
                
            end
        
        end
        
        
        function outputArg = finishUpdate(obj)
            
            obj.x = obj.newStateVariables.xNew;
            obj.y = obj.newStateVariables.yNew;
            obj.u = obj.newStateVariables.uNew;
            obj.w = obj.newStateVariables.wNew;
            
        end
        
        
        
        %% check the positive definiteness of the network version of the obj.testMatrix
        function isPositiveDefinite = checkPositiveDefiniteness(obj, previousSubsystems, subsystems)
            i = length(previousSubsystems)+1;
            iInd = obj.index;
%             disp(['Checking Positive Definiteness at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            
            % We need to make sure W > 0, which happens iff [W_ii, W_i; W_i', M_i] for all i
            W_ii = obj.testMatrix{iInd};
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii
                isPositiveDefinite = all(eig(W_ii)>0);
                tildeW_i = W_ii; % Note that here, \tilde{W}_ii = W_ii = \tilde{W}_i.
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                
                if ~isPositiveDefinite
%                     disp(['Test Matrix is not positive definite at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            else
                % This subsystem has to talk with all the previosSubsystems
                % tildeW_ii>0 iff [W_ii, W_i; W_i', M_i] > 0 is required where 
                % M_i = inv((scriptD_i*scriptA_i)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i)^{-1}')
                
                % W_ii term: known already 
                
                % M_i term
                blockSize = size(W_ii,1); 
                scriptA_i = [];
                scriptD_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % Getting stored info from jInd to create \mathcal{A}_i and \mathcal{D}_i matrices (their j-th columns)
                    tildeW_j = subsystems(jInd).dataToBeDistributed.tildeW;
                                        
                    Z = zeros(blockSize*(i-1-j),blockSize); % (i-1)-j blocks of blockSizeXblockSize zero matrices
                    z = zeros(blockSize*(j-1),blockSize);
                    if j==1
                        tildeW_jj = tildeW_j;                    
                        scriptA_i = [tildeW_jj; Z];         % The first column of \mathcal{A}_i.
                        scriptD_i = [inv(tildeW_jj); Z];    % The first column of \mathcal{D}_i.
                    else
                        tildeW_jj = tildeW_j(:,blockSize*(j-1)+1:blockSize*j);   % last blockSizeXblockSize block in the row block vector
                        tildeW_j  = tildeW_j(:,1:blockSize*(j-1));                % first (j-1) blocks in the row block vector
                        scriptA_i = [scriptA_i, [tildeW_j'; tildeW_jj ; Z]];    % The j-th column of \mathcal{A}_i.
                        scriptD_i = [scriptD_i, [z; inv(tildeW_jj); Z]];         % The j-th column of \mathcal{D}_i.
                    end                    
                end
%                 disp(['Data at ',num2str(iInd),' after ',num2str(previousSubsystems),'.'])
%                 scriptD_i
%                 scriptA_i
                M1_i = inv(scriptD_i*scriptA_i);
                M_i = scriptA_i'*scriptD_i'*scriptA_i;
                if issymmetric(scriptD_i) & issymmetric(scriptA_i) & ~issymmetric(M_i)
                    tf = norm(M_i-M_i.',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M_i = 0.5*(M_i + M_i');
                end
                
                % W_i term
                W_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % W_ij term
                    W_ij = obj.testMatrix{jInd};
                    W_i = [W_i, W_ij];
                end
                
                % Testing [W_ii, W_i; W_i', M_i] > 0
                testMat = [W_ii, W_i; W_i', M_i];
                isPositiveDefinite = all(eig(testMat)>0);
                                
                tildeW_i = W_i*M1_i;
                tildeW_ii = W_ii - tildeW_i*scriptD_i*tildeW_i';
                % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                tildeW_i = [tildeW_i, tildeW_ii];

                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                
                if ~isPositiveDefinite
%                     disp(['Test Matrix is not positive definite at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            
            end
        end
        
        
        
        %% Local stability analysis
        
        function isFeasible = stabilityAnalysis(obj,previousSubsystems,subsystems,solverOptions)
            
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Checking stability at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            A_ii = obj.A{iInd};
            n_i = size(A_ii,1);
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii
                % W_ii = -A_ii^T P_i - P_i A_ii
                % we need M = W_ii>0 and P_i > 0
                
                P_ii = sdpvar(n_i,n_i);
                W_ii = -A_ii'*P_ii - P_ii*A_ii;
                con1 = P_ii >= 0;
                con2 = W_ii >= 0;
                sol = optimize([con1,con2],[],solverOptions);
                isFeasible = sol.problem==0;
                
                P_iiVal = value(P_ii);
                W_iiVal = value(W_ii);
                tildeW_i = W_iiVal; % Note that, here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                
                disp(['Data saved ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_ii; % Storing
                if ~isFeasible
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            
            else
                % This subsystem has to talk with all the previosSubsystems
                % tildeW_ii > 0 iff [M_i, W_i'; W_i, W_ii] > 0 is required where 
                % M_i = inv((scriptD_i*scriptA_i^T)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i^T)^{-1}')
                
                % M_i term
                blockSize = obj.dim_n; 
                scriptA_i = [];
                scriptD_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % Getting stored info from jInd to create \mathcal{A}_i and \mathcal{D}_i matrices (their j-th columns)
                    tildeW_j = subsystems(jInd).dataToBeDistributed.tildeW;
                                        
                    Z = zeros(blockSize, blockSize*(i-1-j)); % (i-1)-j blocks of blockSize X blockSize zero matrices
                    z = zeros(blockSize, blockSize*(j-1));
                    if j==1
                        tildeW_jj = tildeW_j;                    
                        scriptA_i = [tildeW_jj, Z];         % The first row of \mathcal{A}_i.
                        scriptD_i = [inv(tildeW_jj), Z];    % The first row of \mathcal{D}_i.
                    else
                        tildeW_jj = tildeW_j(:,blockSize*(j-1)+1:blockSize*j);   % last blockSizeXblockSize block in the row block vector
                        tildeW_j  = tildeW_j(:,1:blockSize*(j-1));                % first (j-1) blocks in the row block vector
                        scriptA_i = [scriptA_i; [tildeW_j, tildeW_jj, Z]];    % The j-th row of \mathcal{A}_i.
                        scriptD_i = [scriptD_i; [z, inv(tildeW_jj), Z]];         % The j-th row of \mathcal{D}_i.
                    end                    
                end
                disp(['Data at ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                scriptA_i;
                scriptD_i;                
                
                M1_i = inv(scriptD_i*scriptA_i');
                % M_i = inv(M1_i*scriptD_i*M1_i') % THis fills (i-1)x(i-1) blocks in the LMI
                M_i = scriptA_i*scriptD_i*scriptA_i';
                              
                if issymmetric(scriptD_i) & issymmetric(scriptA_i) & ~issymmetric(M_i)
                    tf = norm(M_i-M_i.',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M_i = 0.5*(M_i + M_i');
                end
                
                
                P_ii = sdpvar(n_i,n_i);
                % W_ii and W_i terms
                W_ii = -A_ii'*P_ii - P_ii*A_ii;                
                W_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % W_ij term
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};                    
                    P_jj = subsystems(jInd).dataToBeDistributed.P;
                    n_j = 1+(j-1)+1;
                    
                    W_ij = -A_ji'*P_jj - P_ii*A_ij;
                    W_i = [W_i, W_ij];
                end
                
                con1 = P_ii >= 0;
                con2 = [M_i, W_i';W_i, W_ii] >= 0;
                sol = optimize([con1,con2],[],solverOptions);
                isFeasible = sol.problem==0;
                P_iiVal = value(P_ii); % This needs to be stored
                W_iVal = value(W_i);
                W_iiVal = value(W_ii);
                
                % Need to compute \tilede{W}_i and \tilde{W}_{ii} for storage
                tildeW_i = W_iVal*M1_i;
                tildeW_ii = W_iiVal - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                
                disp(['Data saved ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                tildeW_i = [tildeW_i, tildeW_ii];

                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_ii; % Storing
                if ~isFeasible
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
                
            end
            
        end
        
        
        
        
        
        %% local FSF stabilization
        
        function [isStabilizable,K_ii,K_ijVals,K_jiVals] = FSFStabilization(obj,previousSubsystems, subsystems, solverOptions)
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Stabilizing at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            % W_ij = - M_ii A_ji^T - A_ij M_jj - L_ji^T B_jj^T - B_ii L_ij
            % L_ij is p_i x n_j  and  K_ij = L_ij M_jj^{-1}
            
            A_ii = obj.A{iInd};
            B_ii = obj.B{iInd};
            n_i = obj.dim_n;
            p_i = obj.dim_p; 
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii > 0  
                
                M_ii = sdpvar(n_i,n_i);
                L_ii = sdpvar(p_i,n_i,'full');
                
                W_ii = - M_ii*A_ii' - A_ii*M_ii - L_ii'*B_ii' - B_ii*L_ii;
                con1 = M_ii >= 0;
                con2 = W_ii >= 0;
                sol = optimize([con1,con2],[],solverOptions);
                isStabilizable = sol.problem==0;
                
                M_iiVal = value(M_ii);
                L_iiVal = value(L_ii);
                W_iiVal = value(W_ii);
                tildeW_i = W_iiVal; % Note that, here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                
                
                K_ii = L_iiVal/M_iiVal; % This needs to be stored
                K_jiVals = [];
                K_ijVals = [];
                
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.M = M_iiVal; % Storing
                obj.controllerGains.decenFSFStabCont{iInd} = K_ii;
                
                disp(['Data saved ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                if ~isStabilizable
                    disp(['Not stabilizable at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            
            else
                % This subsystem has to talk with all the previosSubsystems
                % tildeW_ii > 0 iff [M_i, W_i'; W_i, W_ii] > 0 is required where 
                % M_i =
                % inv((scriptD_i*scriptA_i^T)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i^T)^{-1}') = scriptA_i*scriptD_i*scriptA_i'
                
                
                % M_i term
                blockSize = obj.dim_n; 
                scriptA_i = [];
                scriptD_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % Getting stored info from jInd to create \mathcal{A}_i and \mathcal{D}_i matrices (their j-th columns)
                    tildeW_j = subsystems(jInd).dataToBeDistributed.tildeW;
                                        
                    Z = zeros(blockSize, blockSize*(i-1-j)); % (i-1)-j blocks of blockSize X blockSize zero matrices
                    z = zeros(blockSize, blockSize*(j-1));
                    if j==1
                        tildeW_jj = tildeW_j;                    
                        scriptA_i = [tildeW_jj, Z];         % The first row of \mathcal{A}_i.
                        scriptD_i = [inv(tildeW_jj), Z];    % The first row of \mathcal{D}_i.
                    else
                        tildeW_jj = tildeW_j(:,blockSize*(j-1)+1:blockSize*j);   % last blockSizeXblockSize block in the row block vector
                        tildeW_j  = tildeW_j(:,1:blockSize*(j-1));                % first (j-1) blocks in the row block vector
                        scriptA_i = [scriptA_i; [tildeW_j, tildeW_jj, Z]];    % The j-th row of \mathcal{A}_i.
                        scriptD_i = [scriptD_i; [z, inv(tildeW_jj), Z]];         % The j-th row of \mathcal{D}_i.
                    end                    
                end
                disp(['Data at ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                scriptA_i;
                scriptD_i;                
                
                M1_i = inv(scriptD_i*scriptA_i');
                % M_i = inv(M1_i*scriptD_i*M1_i') % THis fills (i-1)x(i-1) blocks in the LMI
                M_i = scriptA_i*scriptD_i*scriptA_i';
                              
                if issymmetric(scriptD_i) & issymmetric(scriptA_i) & ~issymmetric(M_i)
                    tf = norm(M_i-M_i.',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M_i = 0.5*(M_i + M_i');
                end
                
                
                % W_ii and W_i terms
                M_ii = sdpvar(n_i,n_i);
                L_ii = sdpvar(p_i,n_i,'full');
                W_ii = - M_ii*A_ii' - A_ii*M_ii - L_ii'*B_ii' - B_ii*L_ii;
                W_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    n_j = subsystems(jInd).dim_n;
                    p_j = subsystems(jInd).dim_p;
                    
                    if any(subsystems(iInd).neighbors==jInd)
                        L_ij{j} = sdpvar(p_i,n_j,'full');
                    else
                        L_ij{j} = zeros(p_i,n_j);
                    end
                    if any(subsystems(jInd).neighbors==iInd) 
                        L_ji{j} = sdpvar(p_j,n_i,'full');
                    else
                        L_ji{j} = zeros(p_j,n_i);
                    end
                    
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};
                    B_jj  = subsystems(jInd).B{jInd};
                    M_jj = subsystems(jInd).dataToBeDistributed.M;
                    
                    % W_ij = - M_ii A_ji^T - A_ij M_jj - L_ji^T B_jj^T - B_ii L_ij
                    W_ij = - M_ii*A_ji' - A_ij*M_jj - L_ji{j}'*B_jj' - B_ii*L_ij{j};
                    W_i = [W_i, W_ij];
                end
                   
                con1 = M_ii >= 0;
                con2 = [M_i, W_i';W_i, W_ii] >= 0;
                sol = optimize([con1,con2],[],solverOptions);
                isStabilizable = sol.problem==0;
                M_iiVal = value(M_ii); % This needs to be stored
                L_iiVal = value(L_ii);
                W_iVal = value(W_i);
                W_iiVal = value(W_ii);
                
                K_ii = L_iiVal/M_iiVal;
                obj.controllerGains.decenFSFStabCont{iInd} = K_ii;
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    M_jj = subsystems(jInd).dataToBeDistributed.M;
                    
                    K_ij = value(L_ij{j})/M_jj;
                    K_ji = value(L_ji{j})/M_iiVal;
                    K_ijVals{jInd} = K_ij;
                    obj.controllerGains.decenFSFStabCont{jInd} = K_ij;
                    K_jiVals{jInd} = K_ji; % these values will be loaded outside the function
                end  
                
                % Need to compute \tilede{W}_i and \tilde{W}_{ii} for storage
                tildeW_i = W_iVal*M1_i;
                tildeW_ii = W_iiVal - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                
                disp(['Data saved ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                tildeW_i = [tildeW_i, tildeW_ii];
                
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.M = M_iiVal; % Storing
                if ~isStabilizable
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
                
            end
        end
        
        
        
        %% local stable observer design
        
        function [isObserverStable,L_ii,L_ijVals,L_jiVals] = stableObserverDesign(obj,previousSubsystems, subsystems, solverOptions)
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Stabilizing at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            % W_ij = - A_ji^T P_jj - P_ii A_ij + C_ii^T K_ji^T + K_ij C_jj
            % K_ij is n_i x m_j  and  L_ij = P_ii^{-1}K_ij
            
            A_ii = obj.A{iInd};
            C_ii = obj.C{iInd};
            n_i = obj.dim_n;
            m_i = obj.dim_m; 
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii > 0  
                
                P_ii = sdpvar(n_i,n_i);
                K_ii = sdpvar(n_i,m_i,'full');
                
                W_ii = -A_ii'*P_ii - P_ii*A_ii + C_ii'*K_ii' + K_ii*C_ii;
                con1 = P_ii >= 0;
                con2 = W_ii >= 0;
                sol = optimize([con1,con2],[],solverOptions);
                isObserverStable = sol.problem==0;
                
                P_iiVal = value(P_ii);
                K_iiVal = value(K_ii);
                W_iiVal = value(W_ii);
                tildeW_i = W_iiVal; % Note that, here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                
                
                L_ii = P_iiVal\K_iiVal; % This needs to be stored
                L_jiVals = [];
                L_ijVals = [];
                
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_iiVal; % Storing
                obj.observerGains.decenStabObs{iInd} = L_ii;
                
                disp(['Data saved ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                if ~isObserverStable
                    disp(['Not stabilizable at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            
            else
                % This subsystem has to talk with all the previosSubsystems
                % tildeW_ii > 0 iff [M_i, W_i'; W_i, W_ii] > 0 is required where 
                % inv((scriptD_i*scriptA_i^T)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i^T)^{-1}') = scriptA_i*scriptD_i*scriptA_i'
                
                
                % M_i term
                blockSize = obj.dim_n; 
                scriptA_i = [];
                scriptD_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % Getting stored info from jInd to create \mathcal{A}_i and \mathcal{D}_i matrices (their j-th columns)
                    tildeW_j = subsystems(jInd).dataToBeDistributed.tildeW;
                                        
                    Z = zeros(blockSize, blockSize*(i-1-j)); % (i-1)-j blocks of blockSize X blockSize zero matrices
                    z = zeros(blockSize, blockSize*(j-1));
                    if j==1
                        tildeW_jj = tildeW_j;                    
                        scriptA_i = [tildeW_jj, Z];         % The first row of \mathcal{A}_i.
                        scriptD_i = [inv(tildeW_jj), Z];    % The first row of \mathcal{D}_i.
                    else
                        tildeW_jj = tildeW_j(:,blockSize*(j-1)+1:blockSize*j);   % last blockSizeXblockSize block in the row block vector
                        tildeW_j  = tildeW_j(:,1:blockSize*(j-1));                % first (j-1) blocks in the row block vector
                        scriptA_i = [scriptA_i; [tildeW_j, tildeW_jj, Z]];    % The j-th row of \mathcal{A}_i.
                        scriptD_i = [scriptD_i; [z, inv(tildeW_jj), Z]];         % The j-th row of \mathcal{D}_i.
                    end                    
                end
                disp(['Data at ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                scriptA_i;
                scriptD_i;                
                
                M1_i = inv(scriptD_i*scriptA_i');
                % M_i = inv(M1_i*scriptD_i*M1_i') % THis fills (i-1)x(i-1) blocks in the LMI
                M_i = scriptA_i*scriptD_i*scriptA_i';
                              
                if issymmetric(scriptD_i) & issymmetric(scriptA_i) & ~issymmetric(M_i)
                    tf = norm(M_i-M_i.',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M_i = 0.5*(M_i + M_i');
                end
                
                
                % W_ii and W_i terms
                P_ii = sdpvar(n_i,n_i);
                K_ii = sdpvar(n_i,m_i,'full');
                W_ii = - A_ii'*P_ii - P_ii*A_ii + C_ii'*K_ii' + K_ii*C_ii;;
                W_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    n_j = subsystems(jInd).dim_n;
                    m_j = subsystems(jInd).dim_m;
                    
                    if any(subsystems(iInd).neighbors==jInd)
                        K_ij{j} = sdpvar(n_i,m_j,'full');
                    else
                        K_ij{j} = zeros(n_i,m_j);
                    end
                    if any(subsystems(jInd).neighbors==iInd) 
                        K_ji{j} = sdpvar(n_j,m_i,'full');
                    else
                        K_ji{j} = zeros(n_j,m_i);
                    end
                    
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};
                    C_jj = subsystems(jInd).C{jInd};
                    P_jj = subsystems(jInd).dataToBeDistributed.P;
                    
                    % W_ij = - A_ji^T P_jj - P_ii A_ij + C_ii^T K_ji^T + K_ij C_jj
                    W_ij = - A_ji'*P_jj - P_ii*A_ij + C_ii'*K_ji{j}' + K_ij{j}*C_jj;
                    W_i = [W_i, W_ij];
                end
                   
                con1 = P_ii >= 0;
                con2 = [M_i, W_i';W_i, W_ii] >= 0;
                sol = optimize([con1,con2],[],solverOptions);
                isObserverStable = sol.problem==0;
                P_iiVal = value(P_ii); % This needs to be stored
                K_iiVal = value(K_ii);
                W_iVal = value(W_i);
                W_iiVal = value(W_ii);
                
                L_ii = P_iiVal\K_iiVal;
                obj.observerGains.decenStabObs{iInd} = L_ii;
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    P_jj = subsystems(jInd).dataToBeDistributed.P;
                    
                    L_ij = P_iiVal\value(K_ij{j});
                    L_ji = P_jj\value(K_ji{j});
                    L_ijVals{jInd} = L_ij;
                    obj.observerGains.decenStabObs{jInd} = L_ij;
                    L_jiVals{jInd} = L_ji; % these values will be loaded outside the function
                end  
                
                % Need to compute \tilede{W}_i and \tilde{W}_{ii} for storage
                tildeW_i = W_iVal*M1_i;
                tildeW_ii = W_iiVal - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                
                disp(['Data saved ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                tildeW_i = [tildeW_i, tildeW_ii];
                
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_iiVal; % Storing
                if ~isObserverStable
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
                
            end
        end
        
        
        
        %% local DOF stabilization
        
        function [isStabilizable,AcVals,BcVals,CcVals,DcVals] = DOFStabilization(obj,previousSubsystems, subsystems, solverOptions)
            % Note that: AcVals{1} = Ac_ii; AcVals{2} = Ac_ijVals, AcVals{3} = Ac_jiVals
            
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Stabilizing at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            % W1_ij = [Y_ii e_ij,  e_ij;  e_ij,  X_ii e_ij]
            % W2_ij = [W2_11_ij,  W2_12_ij;  W2_21_ij,  W2_22_ij] where
            % W2_11_ij = - A_ij Y_jj - B_ii Cn_ij - Y_ii A_ji^T - Cn_ji^T B_jj^T
            % W2_12_ij = - A_ij - An_ji^T - B_ii Dn_ij C_jj
            % W2_21_ij = - A_ji^T - An_ij - C_ii^T Dn_ji^T B_jj^T
            % W2_22_ij = - X_ii A_ij - Bn_ij C_jj - Aji^T X_jj - C_ii^T Bn_ji^T
            
            % X_ii, Y_ii are n_i x n_i     AND     An_ij is n_i x n_j     AND   Bn_ij is n_i x m_j, 
            % Cn_ij is p_i x n_j           AND     Dn_ij is p_i x m_j 
            
            A_ii = obj.A{iInd};
            B_ii = obj.B{iInd};
            C_ii = obj.C{iInd};
            
            n_i = obj.dim_n;
            m_i = obj.dim_m;
            p_i = obj.dim_p;
            I = eye(n_i,n_i);
            
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W1_ii > 0  and W2_ii > 0
                
                X_ii = sdpvar(n_i,n_i);
                Y_ii = sdpvar(n_i,n_i);
                An_ii = sdpvar(n_i,n_i,'full');
                Bn_ii = sdpvar(n_i,m_i,'full');
                Cn_ii = sdpvar(p_i,n_i,'full');
                Dn_ii = sdpvar(p_i,m_i,'full');
                                
                W1_ii = [Y_ii,  I; I, X_ii];
                W2_11_ii = - A_ii*Y_ii - B_ii*Cn_ii - Y_ii*A_ii' - Cn_ii'*B_ii';
                W2_12_ii = - A_ii - An_ii' - B_ii*Dn_ii*C_ii;
                W2_21_ii = - A_ii' - An_ii - C_ii'*Dn_ii'*B_ii';
                W2_22_ii = - X_ii*A_ii - Bn_ii*C_ii - A_ii'*X_ii - C_ii'*Bn_ii';
                W2_ii = [W2_11_ii, W2_12_ii; W2_21_ii, W2_22_ii];
                
                con1 = X_ii >= 0;
                con2 = Y_ii >= 0;
                con3 = W1_ii >= 0;
                con4 = W2_ii >= 0;
                
                sol = optimize([con1,con2,con3,con4],[],solverOptions);
                isStabilizable = sol.problem==0;
                                
                X_iiVal = value(X_ii);
                Y_iiVal = value(Y_ii);
                An_iiVal = value(An_ii);
                Bn_iiVal = value(Bn_ii);
                Cn_iiVal = value(Cn_ii);
                Dn_iiVal = value(Dn_ii);
                
                W1_iiVal = value(W1_ii);
                W2_iiVal = value(W2_ii);
                tildeW1_i = W1_iiVal; % Note that, here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                tildeW2_i = W2_iiVal; % Note that, here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                
                % CoV:
                % [M_ii,N_ii] = lu(I-X_ii*Y_ii);
                % N_ii = N_ii';
                % Dc_ij = Dn_ij; 
                % Cc_ij = (Cn_ij - Dn_ij C_jj Y_jj)(N_jj')^{-1};
                % Bc_ij = M_ii^{-1}(Bn_ij - X_ii B_ii Dn_ij); 
                % Ac_ij = M_ii^{-1}(An_ij - Bn_ij C_jj Y_jj - X_ii B_ii Cn_ij - X_ii(A_ij - B_ii Dn_ij C_jj)Y_jj)(N_jj')^{-1}; 
                
                [M_ii,N_ii] = lu(I-X_iiVal*Y_iiVal);
                N_ii = N_ii';
                Dc_ii = Dn_iiVal; 
                Cc_ii = (Cn_iiVal - Dn_iiVal*C_ii*Y_iiVal)/(N_ii');
                Bc_ii = M_ii\(Bn_iiVal - X_iiVal*B_ii*Dn_iiVal); 
                Ac_ii = M_ii\(An_iiVal - Bn_iiVal*C_ii*Y_iiVal - X_iiVal*B_ii*Cn_iiVal - X_iiVal*(A_ii - B_ii*Dn_iiVal*C_ii)*Y_iiVal)/(N_ii'); 
                
                obj.controllerGains.decenDOFStabContAc{iInd} = Ac_ii;
                obj.controllerGains.decenDOFStabContBc{iInd} = Bc_ii;
                obj.controllerGains.decenDOFStabContCc{iInd} = Cc_ii;
                obj.controllerGains.decenDOFStabContDc{iInd} = Dc_ii;
                
                % Note that: AcVals{1} = Ac_ii; AcVals{2} = Ac_ijVals, AcVals{3} = Ac_jiVals
                AcVals{1} = Ac_ii; AcVals{2} = []; AcVals{3} = [];
                BcVals{1} = Bc_ii; BcVals{2} = []; BcVals{3} = [];
                CcVals{1} = Cc_ii; CcVals{2} = []; CcVals{3} = [];
                DcVals{1} = Dc_ii; DcVals{2} = []; DcVals{3} = [];
                                
                obj.dataToBeDistributed.X = X_iiVal; % Storing
                obj.dataToBeDistributed.Y = Y_iiVal; % Storing
                obj.dataToBeDistributed.M = M_ii; % Storing
                obj.dataToBeDistributed.N = N_ii; % Storing
                obj.dataToBeDistributed.tildeW1 = tildeW1_i; % Storing
                obj.dataToBeDistributed.tildeW2 = tildeW2_i; % Storing
                
                disp(['Data saved ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                if ~isStabilizable
                    disp(['Not stabilizable at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            
            else
                % This subsystem has to talk with all the previosSubsystems
                % tildeW1_ii > 0 iff [M1_i, W1_i'; W1_i, W1_ii] > 0 is required where 
                % M1_i = scriptA1_i*scriptD1_i*scriptA1_i'; Also: 
                % M11_i = inv(scriptD1_i*scriptA1_i') and tildeW1_i = W1_i*M11_i;
                % Note that: inv((scriptD_i*scriptA_i^T)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i^T)^{-1}') = scriptA_i*scriptD_i*scriptA_i'
                
                
                % M1_i term
                blockSize = 2*n_i; % assume n_i = n_j 
                scriptA1_i = []; scriptA2_i = [];
                scriptD1_i = []; scriptD2_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % Getting stored info from jInd to create \mathcal{A}_i and \mathcal{D}_i matrices (their j-th columns)
                    tildeW1_j = subsystems(jInd).dataToBeDistributed.tildeW1;
                    tildeW2_j = subsystems(jInd).dataToBeDistributed.tildeW2;
                                        
                    Z = zeros(blockSize, blockSize*(i-1-j)); % (i-1)-j blocks of blockSize X blockSize zero matrices
                    z = zeros(blockSize, blockSize*(j-1));
                    if j==1
                        tildeW1_jj = tildeW1_j;
                        tildeW2_jj = tildeW2_j;
                        scriptA1_i = [tildeW1_jj, Z];         % The first row of \mathcal{A}_i.
                        scriptA2_i = [tildeW2_jj, Z];         % The first row of \mathcal{A}_i.
                        scriptD1_i = [inv(tildeW1_jj), Z];    % The first row of \mathcal{D}_i.
                        scriptD2_i = [inv(tildeW2_jj), Z];    % The first row of \mathcal{D}_i.
                    else
                        tildeW1_jj = tildeW1_j(:,blockSize*(j-1)+1:blockSize*j);   % last blockSizeXblockSize block in the row block vector
                        tildeW2_jj = tildeW2_j(:,blockSize*(j-1)+1:blockSize*j);   % last blockSizeXblockSize block in the row block vector
                        tildeW1_j  = tildeW1_j(:,1:blockSize*(j-1));                % first (j-1) blocks in the row block vector
                        tildeW2_j  = tildeW2_j(:,1:blockSize*(j-1));                % first (j-1) blocks in the row block vector
                        scriptA1_i = [scriptA1_i; [tildeW1_j, tildeW1_jj, Z]];    % The j-th row of \mathcal{A}_i.
                        scriptA2_i = [scriptA2_i; [tildeW2_j, tildeW2_jj, Z]];    % The j-th row of \mathcal{A}_i.
                        scriptD1_i = [scriptD1_i; [z, inv(tildeW1_jj), Z]];         % The j-th row of \mathcal{D}_i.
                        scriptD2_i = [scriptD2_i; [z, inv(tildeW2_jj), Z]];         % The j-th row of \mathcal{D}_i.
                    end                    
                end
                disp(['Data at ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                scriptA1_i; scriptA2_i; 
                scriptD1_i; scriptD2_i;                
                
                M11_i = inv(scriptD1_i*scriptA1_i');
                M21_i = inv(scriptD2_i*scriptA2_i');
                M1_i = scriptA1_i*scriptD1_i*scriptA1_i';
                M2_i = scriptA2_i*scriptD2_i*scriptA2_i';
                              
                if issymmetric(scriptD1_i) & issymmetric(scriptA1_i) & ~issymmetric(M1_i)
                    tf = norm(M1_i-M1_i.',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M1_i = 0.5*(M1_i + M1_i');
                end
                
                if issymmetric(scriptD2_i) & issymmetric(scriptA2_i) & ~issymmetric(M2_i)
                    tf = norm(M2_i-M2_i.',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M2_i = 0.5*(M2_i + M2_i');
                end
                
                
                % W1_ii, W1_i and W2_ii, W2_i terms
                X_ii = sdpvar(n_i,n_i);
                Y_ii = sdpvar(n_i,n_i);
                An_ii = sdpvar(n_i,n_i,'full');
                Bn_ii = sdpvar(n_i,m_i,'full');
                Cn_ii = sdpvar(p_i,n_i,'full');
                Dn_ii = sdpvar(p_i,m_i,'full');
                                
                W1_ii = [Y_ii,  I; I, X_ii];
                W2_11_ii = - A_ii*Y_ii - B_ii*Cn_ii - Y_ii*A_ii' - Cn_ii'*B_ii';
                W2_12_ii = - A_ii - An_ii' - B_ii*Dn_ii*C_ii;
                W2_21_ii = - A_ii' - An_ii - C_ii'*Dn_ii'*B_ii';
                W2_22_ii = - X_ii*A_ii - Bn_ii*C_ii - A_ii'*X_ii - C_ii'*Bn_ii';
                W2_ii = [W2_11_ii, W2_12_ii; W2_21_ii, W2_22_ii];
                
                W1_i = [];
                W2_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    % An_ij is n_i x n_j    AND   Bn_ij is n_i x m_j, 
                    % Cn_ij is p_i x n_j    AND   Dn_ij is p_i x m_j 
                    
                    m_j = subsystems(jInd).dim_m;
                    n_j = subsystems(jInd).dim_n;
                    p_j = subsystems(jInd).dim_p;
                    
                    if any(subsystems(iInd).neighbors==jInd)
                        An_ij{j} = sdpvar(n_i,n_j,'full');
                        Bn_ij{j} = sdpvar(n_i,m_j,'full');
                        Cn_ij{j} = sdpvar(p_i,n_j,'full');
                        Dn_ij{j} = sdpvar(p_i,m_j,'full');
                    else
                        An_ij{j} = zeros(n_i,n_j);
                        Bn_ij{j} = zeros(n_i,m_j);
                        Cn_ij{j} = zeros(p_i,n_j);
                        Dn_ij{j} = zeros(p_i,m_j);
                    end
                    if any(subsystems(jInd).neighbors==iInd)
                        An_ji{j} = sdpvar(n_j,n_i,'full');
                        Bn_ji{j} = sdpvar(n_j,m_i,'full');
                        Cn_ji{j} = sdpvar(p_j,n_i,'full');
                        Dn_ji{j} = sdpvar(p_j,m_i,'full');
                    else
                        An_ji{j} = zeros(n_j,n_i);
                        Bn_ji{j} = zeros(n_j,m_i);
                        Cn_ji{j} = zeros(p_j,n_i);
                        Dn_ji{j} = zeros(p_j,m_i);
                    end
                    
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};
                    B_jj = subsystems(jInd).B{jInd};
                    C_jj = subsystems(jInd).C{jInd};
                    X_jj = subsystems(jInd).dataToBeDistributed.X;
                    Y_jj = subsystems(jInd).dataToBeDistributed.Y;
                    
                    
                    % W1_ij = [Y_ii e_ij,  e_ij;  e_ij,  X_ii e_ij]
                    W1_ij = [Y_ii*(i==j),  I*(i==j);  I*(i==j),  X_ii*(i==j)];
                    
                    % W2_11_ij = - A_ij Y_jj - B_ii Cn_ij - Y_ii A_ji^T - Cn_ji^T B_jj^T
                    W2_11_ij = - A_ij*Y_jj - B_ii*Cn_ij{j} - Y_ii*A_ji' - Cn_ji{j}'*B_jj';
                    % W2_12_ij = - A_ij - An_ji^T - B_ii Dn_ij C_jj
                    W2_12_ij = - A_ij - An_ji{j}' - B_ii*Dn_ij{j}*C_jj;
                    % W2_21_ij = - A_ji^T - An_ij - C_ii^T Dn_ji^T B_jj^T
                    W2_21_ij = - A_ji' - An_ij{j} - C_ii'*Dn_ji{j}'*B_jj';
                    % W2_22_ij = - X_ii A_ij - Bn_ij C_jj - Aji^T X_jj - C_ii^T Bn_ji^T
                    W2_22_ij = - X_ii*A_ij - Bn_ij{j}*C_jj - A_ji'*X_jj - C_ii'*Bn_ji{j}';
                    % W2_ij = [W2_11_ij,  W2_12_ij;  W2_21_ij,  W2_22_ij]                     
                    W2_ij = [W2_11_ij,  W2_12_ij;  W2_21_ij,  W2_22_ij];
                    
                    W1_i = [W1_i, W1_ij];
                    W2_i = [W2_i, W2_ij];
                end
                                
                con1 = X_ii >= 0;
                con2 = Y_ii >= 0;
                con3 = [M1_i, W1_i';W1_i, W1_ii] >= 0;
                con4 = [M2_i, W2_i';W2_i, W2_ii] >= 0;
                sol = optimize([con1,con2,con3,con4],[],solverOptions);
                isStabilizable = sol.problem==0;
                
                X_iiVal = value(X_ii);
                Y_iiVal = value(Y_ii);
                An_iiVal = value(An_ii);
                Bn_iiVal = value(Bn_ii);
                Cn_iiVal = value(Cn_ii);
                Dn_iiVal = value(Dn_ii);
                
                W1_iVal = value(W1_i);
                W2_iVal = value(W2_i);
                W1_iiVal = value(W1_ii);
                W2_iiVal = value(W2_ii);
                
                
                % CoV:
                % [M_ii,N_ii] = lu(I-X_ii*Y_ii);
                % N_ii = N_ii';
                % Dc_ij = Dn_ij; 
                % Cc_ij = (Cn_ij - Dn_ij C_jj Y_jj)(N_jj')^{-1};
                % Bc_ij = M_ii^{-1}(Bn_ij - X_ii B_ii Dn_ij); 
                % Ac_ij = M_ii^{-1}(An_ij - Bn_ij C_jj Y_jj - X_ii B_ii Cn_ij - X_ii(A_ij - B_ii Dn_ij C_jj)Y_jj)(N_jj')^{-1}; 
                
                [M_ii,N_ii] = lu(I-X_iiVal*Y_iiVal);
                N_ii = N_ii';
                Dc_ii = Dn_iiVal; 
                Cc_ii = (Cn_iiVal - Dn_iiVal*C_ii*Y_iiVal)/(N_ii');
                Bc_ii = M_ii\(Bn_iiVal - X_iiVal*B_ii*Dn_iiVal); 
                Ac_ii = M_ii\(An_iiVal - Bn_iiVal*C_ii*Y_iiVal - X_iiVal*B_ii*Cn_iiVal - X_iiVal*(A_ii - B_ii*Dn_iiVal*C_ii)*Y_iiVal)/(N_ii'); 
                
                obj.controllerGains.decenDOFStabContAc{iInd} = Ac_ii;
                obj.controllerGains.decenDOFStabContBc{iInd} = Bc_ii;
                obj.controllerGains.decenDOFStabContCc{iInd} = Cc_ii;
                obj.controllerGains.decenDOFStabContDc{iInd} = Dc_ii;
                
                % Note that: AcVals{1} = Ac_ii; AcVals{2} = Ac_ijVals, AcVals{3} = Ac_jiVals
                AcVals{1} = Ac_ii; 
                BcVals{1} = Bc_ii; 
                CcVals{1} = Cc_ii; 
                DcVals{1} = Dc_ii; 
                
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    X_jj = subsystems(jInd).dataToBeDistributed.X;
                    Y_jj = subsystems(jInd).dataToBeDistributed.Y;
                    M_jj = subsystems(jInd).dataToBeDistributed.M;
                    N_jj = subsystems(jInd).dataToBeDistributed.N;
                    
                    % CoV:
                    % Dc_ij = Dn_ij; 
                    Dc_ij = value(Dn_ij{j});
                    Dc_ji = value(Dn_ji{j});
                    % Cc_ij = (Cn_ij - Dn_ij C_jj Y_jj)(N_jj')^{-1};
                    Cc_ij = value((Cn_ij{j} - Dn_ij{j}*C_jj*Y_jj)/(N_jj'));
                    Cc_ji = value((Cn_ji{j} - Dn_ji{j}*C_ii*Y_iiVal)/(N_ii'));
                    % Bc_ij = M_ii^{-1}(Bn_ij - X_ii B_ii Dn_ij); 
                    Bc_ij = value(M_ii\(Bn_ij{j} - X_iiVal*B_ii*Dn_ij{j}));
                    Bc_ji = value(M_jj\(Bn_ji{j} - X_jj*B_jj*Dn_ji{j}));
                    % Ac_ij = M_ii^{-1}(An_ij - Bn_ij C_jj Y_jj - X_ii B_ii Cn_ij - X_ii(A_ij - B_ii Dn_ij C_jj)Y_jj)(N_jj')^{-1}; 
                    Ac_ij = value(M_ii\(An_ij{j} - Bn_ij{j}*C_jj*Y_jj - X_iiVal*B_ii*Cn_ij{j} - X_iiVal*(A_ij - B_ii*Dn_ij{j}*C_jj)*Y_jj)/(N_jj')); 
                    Ac_ji = value(M_jj\(An_ji{j} - Bn_ji{j}*C_ii*Y_iiVal - X_jj*B_jj*Cn_ji{j} - X_jj*(A_ji - B_jj*Dn_ji{j}*C_ii)*Y_iiVal)/(N_ii')); 
                    
                    
                    % Note that: AcVals{1} = Ac_ii; AcVals{2} = Ac_ijVals, AcVals{3} = Ac_jiVals
                    AcVals{2}{jInd} = Ac_ij; AcVals{3}{jInd} = Ac_ji;
                    BcVals{2}{jInd} = Bc_ij; BcVals{3}{jInd} = Bc_ji;
                    CcVals{2}{jInd} = Cc_ij; CcVals{3}{jInd} = Cc_ji;
                    DcVals{2}{jInd} = Dc_ij; DcVals{3}{jInd} = Dc_ji;
                    
                    obj.controllerGains.decenDOFStabContAc{jInd} = Ac_ij;
                    obj.controllerGains.decenDOFStabContBc{jInd} = Bc_ij;
                    obj.controllerGains.decenDOFStabContCc{jInd} = Cc_ij;
                    obj.controllerGains.decenDOFStabContDc{jInd} = Dc_ij;
                    % Ac_ji, Bc_ji, Cc_ji and Dc_ji values will be loaded outside the function
                end  
                
                
                % Need to compute \tilede{W}_i and \tilde{W}_{ii} for storage
                tildeW1_i = W1_iVal*M11_i;
                tildeW2_i = W2_iVal*M21_i;
                tildeW1_ii = W1_iiVal - tildeW1_i*scriptD1_i*tildeW1_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                tildeW2_ii = W2_iiVal - tildeW2_i*scriptD2_i*tildeW2_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                tildeW1_i = [tildeW1_i, tildeW1_ii];
                tildeW2_i = [tildeW2_i, tildeW2_ii];
                
                obj.dataToBeDistributed.tildeW1 = tildeW1_i; % Storing
                obj.dataToBeDistributed.tildeW2 = tildeW2_i; % Storing
                
                obj.dataToBeDistributed.X = X_iiVal; % Storing
                obj.dataToBeDistributed.Y = Y_iiVal; % Storing
                obj.dataToBeDistributed.M = M_ii; % Storing
                obj.dataToBeDistributed.N = N_ii; % Storing
                
                disp(['Data saved ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                if ~isStabilizable
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
                
            end
        end
        
        
        
        
        %% Local QSR-Dissipativity analysis
        function isFeasible = dissipativityAnalysis(obj,dissFrom,dissTo,previousSubsystems,subsystems,solverOptions)
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Checking QSR-Dissipativity at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            
            Q_ii = obj.dataToBeDistributed.Q{iInd};
            S_ii = obj.dataToBeDistributed.S{iInd};
            R_ii = obj.dataToBeDistributed.R{iInd};
            
            A_ii = obj.A{iInd};
            if isequal(dissFrom,'u')  
                B_ii = obj.B{iInd};
                if isequal(dissTo,'y')
                    C_ii = obj.C{iInd};
                    D_ii = obj.D{iInd};
                elseif isequal(dissTo,'z')
                    C_ii = obj.G{iInd};
                    D_ii = obj.H{iInd};
                end                
            elseif isequal(dissFrom,'w')
                B_ii = obj.E{iInd};
                if isequal(dissTo,'y')
                    C_ii = obj.C{iInd};
                    D_ii = obj.F{iInd};
                elseif isequal(dissTo,'z')
                    C_ii = obj.G{iInd};
                    D_ii = obj.J{iInd};
                end
            end
                        
            n_i = size(A_ii,1);
                       
            
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii
                % W_ii = -A_ii^T P_i - P_i A_ii
                % we need M = W_ii>0 and P_i > 0
                
                P_ii = sdpvar(n_i,n_i);
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                W_ii = -A_ii'*P_ii - P_ii*A_ii;
                
                
                
                
                con1 = P_ii >= 0;
                con2 = W_ii >= 0;
                sol = optimize([con1,con2],[],solverOptions);
                isFeasible = sol.problem==0;
                
                P_iiVal = value(P_ii);
                W_iiVal = value(W_ii);
                tildeW_i = W_iiVal; % Note that, here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                
                disp(['Data saved ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_ii; % Storing
                if ~isFeasible
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            
            else
                % This subsystem has to talk with all the previosSubsystems
                % tildeW_ii > 0 iff [M_i, W_i'; W_i, W_ii] > 0 is required where 
                % M_i = inv((scriptD_i*scriptA_i^T)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i^T)^{-1}')
                
                % M_i term
                blockSize = obj.dim_n; 
                scriptA_i = [];
                scriptD_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % Getting stored info from jInd to create \mathcal{A}_i and \mathcal{D}_i matrices (their j-th columns)
                    tildeW_j = subsystems(jInd).dataToBeDistributed.tildeW;
                                        
                    Z = zeros(blockSize, blockSize*(i-1-j)); % (i-1)-j blocks of blockSize X blockSize zero matrices
                    z = zeros(blockSize, blockSize*(j-1));
                    if j==1
                        tildeW_jj = tildeW_j;                    
                        scriptA_i = [tildeW_jj, Z];         % The first row of \mathcal{A}_i.
                        scriptD_i = [inv(tildeW_jj), Z];    % The first row of \mathcal{D}_i.
                    else
                        tildeW_jj = tildeW_j(:,blockSize*(j-1)+1:blockSize*j);   % last blockSizeXblockSize block in the row block vector
                        tildeW_j  = tildeW_j(:,1:blockSize*(j-1));                % first (j-1) blocks in the row block vector
                        scriptA_i = [scriptA_i; [tildeW_j, tildeW_jj, Z]];    % The j-th row of \mathcal{A}_i.
                        scriptD_i = [scriptD_i; [z, inv(tildeW_jj), Z]];         % The j-th row of \mathcal{D}_i.
                    end                    
                end
                disp(['Data at ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                scriptA_i;
                scriptD_i;                
                
                M1_i = inv(scriptD_i*scriptA_i');
                % M_i = inv(M1_i*scriptD_i*M1_i') % THis fills (i-1)x(i-1) blocks in the LMI
                M_i = scriptA_i*scriptD_i*scriptA_i';
                              
                if issymmetric(scriptD_i) & issymmetric(scriptA_i) & ~issymmetric(M_i)
                    tf = norm(M_i-M_i.',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M_i = 0.5*(M_i + M_i');
                end
                
                
                P_ii = sdpvar(n_i,n_i);
                % W_ii and W_i terms
                W_ii = -A_ii'*P_ii - P_ii*A_ii;                
                W_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % W_ij term
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};                    
                    P_jj = subsystems(jInd).dataToBeDistributed.P;
                    n_j = 1+(j-1)+1;
                    
                    W_ij = -A_ji'*P_jj - P_ii*A_ij;
                    W_i = [W_i, W_ij];
                end
                
                con1 = P_ii >= 0;
                con2 = [M_i, W_i';W_i, W_ii] >= 0;
                sol = optimize([con1,con2],[],solverOptions);
                isFeasible = sol.problem==0;
                P_iiVal = value(P_ii); % This needs to be stored
                W_iVal = value(W_i);
                W_iiVal = value(W_ii);
                
                % Need to compute \tilede{W}_i and \tilde{W}_{ii} for storage
                tildeW_i = W_iVal*M1_i;
                tildeW_ii = W_iiVal - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                
                disp(['Data saved ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                tildeW_i = [tildeW_i, tildeW_ii];

                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_ii; % Storing
                if ~isFeasible
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
                
            end
          
        end
        
        
        
        
        
        
        
        
        %% Local dissipative controller design using LMIs
        function [isDissipative,K_ii,K_ijVals,K_jiVals] = designLocalDissipatingSFBControllers(obj,previousSubsystems, subsystems)
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Dissipating at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            
            A_ii = obj.A{iInd};
            B_i = obj.B{iInd};
            C_i = obj.C{iInd};
            E_ii = obj.E{iInd};
            F_i = obj.F{iInd};
            Q_ii = obj.dataToBeDistributed.Q{iInd};
            S_ii = obj.dataToBeDistributed.Q{iInd};
            R_ii = obj.dataToBeDistributed.Q{iInd};
            
            dim_n = obj.dim_n;
            dim_p = obj.dim_p;
            dim_q = obj.dim_q;
            
            
            % For Dissipativity: W_ij = [W_ij1, W_ij2; W_ij3, W_ij4]
            % W_ij1 = -A_{ji}^T P_j -K_{ji}^T B_j^T P_j - P_iA_{ij} - P_i B_i K_{ij} + C_i^T Q_{ij} C_j
            % W_ij2 = -P_i E_{ij} + C_i^T S_{ij} + C_i^T Q_{ij} F_j
            % W_ij3 = -E_{ji}^T P_j + S_{ji}^T C_j + F_i^T Q_{ji}^T C_j  
            % W_ij4 = F_i^T Q_{ij} F_j + (F_i^T S_{ij} + S_{ji}^T F_j) + R_{ij}
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii > 0 and P_i > 0 with K_ii free
                % W_ii = [W_ii1, W_ii2; W_ii3, W_ii4]
                % W_ii1 = -A_{ii}^T P_i -K_{ii}^T B_i^T P_i - P_iA_{ii} - P_i B_i K_{ii} + C_i^T Q_{ii} C_i
                % W_ii2 = -P_i E_{ii} + C_i^T S_{ii} + C_i^T Q_{ii} F_i
                % W_ii3 = -E_{ii}^T P_i + S_{ii}^T C_i + F_i^T Q_{ii}^T C_i  
                % W_ii4 = F_i^T Q_{ii} F_i + (F_i^T S_{ii} + S_{ii}^T F_i) + R_{ii}
            
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                P = lmivar(1,[dim_n, 1]);   % P is a variable, 1: square symmetric, size(A,1) gives the size and 1 gives that Q is a full matrix 
                L = lmivar(2,[dim_n,dim_n])      % L_i = P_i B_i K_ii is a variable, rectangular matrix
                
                % W_ii1 = -A_{ii}^T P_i -K_{ii}^T B_i^T P_i - P_iA_{ii} - P_i B_i K_{ii} + C_i^T Q_{ii} C_i
                % W_ii1 = -A_{ii}^T P_i -L_i^T - P_iA_{ii} - L_i + C_i^T Q_{ii} C_i
                lmiterm([-1, 1, 1, P],-1,A_ii,'s');
                lmiterm([-1, 1, 1, L],-1,1,'s');
                matTemp1 = C_i'*Q_ii*C_i;
                lmiterm([-1, 1, 1, 0],matTemp1);
                % W_ii2 = -P_i E_{ii} + C_i^T S_{ii} + C_i^T Q_{ii} F_i
                lmiterm([-1, 1, 2, P],-1,E_ii); 
                matTemp2 = C_i'*S_ii + C_i'*Q_ii*F_i;
                lmiterm([-1, 1, 2, 0],matTemp2); 
                % W_ii3 = -E_{ii}^T P_i + S_{ii}^T C_i + F_i^T Q_{ii}^T C_i  
                lmiterm([-1, 2, 1, P],-E_ii',1);
                lmiterm([-1, 2, 1, 0],matTemp2'); 
                % W_ii4 = F_i^T Q_{ii} F_i + (F_i^T S_{ii} + S_{ii}^T F_i) + R_{ii}
                matTemp4 = F_i'*Q_ii*F_i + (F_i'*S_ii+S_ii'*F_i) + R_ii;
                lmiterm([-1, 2, 2, 0],matTemp4);
                
                % P>0
                lmiterm([-2, 1, 1, P],1,1); % defines -P<0
                lmisys = getlmis;
                [tmin,sol] = feasp(lmisys); % Solve the LMI system
                isDissipative = tmin <= 0; % strictly feasible if this is satisfied
                
                P_i = dec2mat(lmisys,sol,P); % This needs to be stored
                L_i = dec2mat(lmisys,sol,L); % This needs to be stored
                K_ii = pinv(B_i)*inv(P_i)*L_i;
                K_jiVals = [];
                K_ijVals = [];
                tildeW_ii1 = -A_ii'*P_i - P_i*A_ii + L_i + L_i' + matTemp1;
                tildeW_ii2 = -P_i*E_ii + matTemp2;
                tildeW_ii3 = tildeW_ii2';
                tildeW_ii4 = matTemp4;
                tildeW_i = [tildeW_ii1, tildeW_ii2; tildeW_ii3, tildeW_ii4]; % Note that here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                % LMI End
                
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                obj.controllerGains.localSFBDissLMI{iInd} = K_ii;
                if ~isDissipative
                    disp(['Not dissipative at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
                   
            else
                % This subsystem has to talk with all the previosSubsystems
                % tildeW_ii > 0 iff [W_ii, W_i; W_i', M_i] > 0 is required where 
                % M_i = inv((scriptD_i*scriptA_i)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i)^{-1}')
                
                % Introduce L_i = P_i B_i K_ii;  L_ij = P_i B_i K_ij; L_ji = K_ji' as variables apart from variables P_i.
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                P = lmivar(1,[dim_n, 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    % L_ij = lmivar(2,[dim_n, dim_n])     % For L_ij = P_i B_i K_ij
                    % L_ji = lmivar(2,[dim_n, dim_p])     % For L_ji = K_ji^T
                    eval(['L_i',num2str(j),'= lmivar(2,[dim_n, dim_n]);'])
                    eval(['L_',num2str(j),'i= lmivar(2,[dim_n, dim_p]);'])
                end
                L = lmivar(2,[dim_n,dim_n]); % For L_i = P_i B_i K_ii 
                                
                % W_ii term: W_ii = [W_ii1, W_ii2; W_ii3, W_ii4]
                % W_ii1 = -A_{ii}^T P_i -K_{ii}^T B_i^T P_i - P_iA_{ii} - P_i B_i K_{ii} + C_i^T Q_{ii} C_i
                % W_ii1 = -A_{ii}^T P_i -L_i^T - P_iA_{ii} - L_i + C_i^T Q_{ii} C_i
                lmiterm([-1, 1, 1, P],-1,A_ii,'s');
                lmiterm([-1, 1, 1, L],-1,1,'s');
                matTemp1 = C_i'*Q_ii*C_i;
                lmiterm([-1, 1, 1, 0],matTemp1);
                % W_ii2 = -P_i E_{ii} + C_i^T S_{ii} + C_i^T Q_{ii} F_i
                lmiterm([-1, 1, 2, P],-1,E_ii); 
                matTemp2 = C_i'*S_ii + C_i'*Q_ii*F_i;
                lmiterm([-1, 1, 2, 0],matTemp2); 
                % W_ii3 = -E_{ii}^T P_i + S_{ii}^T C_i + F_i^T Q_{ii}^T C_i  
                lmiterm([-1, 2, 1, P],-E_ii',1);
                lmiterm([-1, 2, 1, 0],matTemp2'); 
                % W_ii4 = F_i^T Q_{ii} F_i + (F_i^T S_{ii} + S_{ii}^T F_i) + R_{ii}
                matTemp4 = F_i'*Q_ii*F_i + (F_i'*S_ii+S_ii'*F_i) + R_ii;
                lmiterm([-1, 2, 2, 0],matTemp4);
                
                
                % M_i term = inv((scriptD_i*scriptA_i)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i)^{-1}')
                blockSize = dim_n+dim_q;
                scriptA_i = [];
                scriptD_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % Getting stored info from jInd to create \mathcal{A}_i and \mathcal{D}_i matrices (their j-th columns)
                    tildeW_j = subsystems(jInd).dataToBeDistributed.tildeW;
                                        
                    Z = zeros(blockSize*(i-1-j),blockSize); % (i-1)-j blocks of blockSizeXblockSize zero matrices
                    z = zeros(blockSize*(j-1),blockSize);
                    if j==1
                        tildeW_jj = tildeW_j;                    
                        scriptA_i = [tildeW_jj; Z];         % The first column of \mathcal{A}_i.
                        scriptD_i = [inv(tildeW_jj); Z];    % The first column of \mathcal{D}_i.
                    else
                        tildeW_jj = tildeW_j(:,blockSize*(j-1)+1:blockSize*j);   % last blockSizeXblockSize block in the row block vector
                        tildeW_j  = tildeW_j(:,1:blockSize*(j-1));                % first (j-1) blocks in the row block vector
                        scriptA_i = [scriptA_i, [tildeW_j'; tildeW_jj ; Z]];    % The j-th column of \mathcal{A}_i.
                        scriptD_i = [scriptD_i, [z; inv(tildeW_jj); Z]];         % The j-th column of \mathcal{D}_i.
                    end                    
                end
                M1_i = inv(scriptD_i*scriptA_i);
                % M_i = inv(M1_i*scriptD_i*M1_i'); % THis fills 2*(i-1)x2*(i-1) blocks in the LMI
                M_i = scriptA_i'*scriptD_i'*scriptA_i;
                if issymmetric(scriptD_i) & issymmetric(scriptA_i) & ~issymmetric(M_i)
                    tf = norm(M_i-M_i.',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M_i = 0.5*(M_i + M_i');
                end
                
                % Filling M_i
                for x=1:1:length(previousSubsystems)
                    m_x = (x-1)*blockSize; % exact index
                    n_x = 2+2*(x-1)+1;           % block index
                    for y=1:1:length(previousSubsystems)
                        m_y = (y-1)*blockSize;  % exact index
                        n_y = 2+2*(y-1)+1;           % block index
                        M_i_xy = M_i(m_x+1:m_x+blockSize, m_y+1:m_y+blockSize);
                        M_i_xy1 = M_i_xy(1:obj.dim_n,1:obj.dim_n);
                        M_i_xy2 = M_i_xy(1:obj.dim_n,obj.dim_n+1:blockSize);
                        M_i_xy3 = M_i_xy(obj.dim_n+1:blockSize,1:obj.dim_n);
                        M_i_xy4 = M_i_xy(obj.dim_n+1:blockSize,obj.dim_n+1:blockSize);
                        lmiterm([-1, n_x, n_y, 0],M_i_xy1);
                        lmiterm([-1, n_x, n_y+1, 0],M_i_xy2);
                        lmiterm([-1, n_x+1, n_y, 0],M_i_xy3);
                        lmiterm([-1, n_x+1, n_y+1, 0],M_i_xy4);
                    end
                end
               
                % W_i term
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % W_ij term
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};
                    B_j = subsystems(jInd).B{jInd};
                    C_j = subsystems(jInd).C{jInd};
                    E_ij = obj.E{jInd};
                    E_ji = subsystems(jInd).E{iInd};
                    F_j = subsystems(jInd).F{jInd};
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    Q_ij = obj.dataToBeDistributed.Q{jInd};
                    Q_ji = subsystems(jInd).dataToBeDistributed.Q{iInd};
                    S_ij = obj.dataToBeDistributed.S{jInd};
                    S_ji = subsystems(jInd).dataToBeDistributed.S{iInd};
                    R_ij = obj.dataToBeDistributed.R{jInd};
                    n_j = 2+2*(j-1)+1;
                    
                    
                    % Filling W_ij = [W_ij1, W_ij2; W_ij3, W_ij4]
                    % W_ij1 = -A_{ji}^T P_j -K_{ji}^T B_j^T P_j - P_iA_{ij} - P_i B_i K_{ij} + C_i^T Q_{ij} C_j
                    % W_ij1 = -A_{ji}^T P_j -L_{ji} B_j^T P_j - P_iA_{ij} - L_{ij} + C_i^T Q_{ij} C_j
                    lmiterm([-1, 1, n_j, P],-1,A_ij);
                    matTemp1 = -A_ji'*P_j + C_i'*Q_ij*C_j;
                    lmiterm([-1, 1, n_j, 0],matTemp1);
                    matTemp2 = B_j'*P_j;
                    % lmiterm([-1, 1, n_j, L_ji],-1,matTemp2);
                    % lmiterm([-1, 1, n_j, L_ij],-1,1);
                    eval(['lmiterm([-1, 1, n_j, L_',num2str(j),'i],-1,matTemp2);']);
                    eval(['lmiterm([-1, 1, n_j, L_i',num2str(j),'],-1,1);']);
                    % W_ij2 = -P_i E_{ij} + C_i^T S_{ij} + C_i^T Q_{ij} F_j
                    lmiterm([-1, 1, n_j+1, P],-1,E_ij); 
                    matTemp3 = C_i'*S_ij + C_i'*Q_ij*F_j;
                    lmiterm([-1, 1, n_j+1, 0],matTemp3);
                    % W_ij3 = -E_{ji}^T P_j + S_{ji}^T C_j + F_i^T Q_{ji}^T C_j  
                    matTemp4 = -E_ji'*P_j + S_ji'*C_j + F_i'*Q_ji'*C_j;
                    lmiterm([-1, 2, n_j, 0],matTemp4);                    
                    % W_ij4 = F_i^T Q_{ij} F_j + (F_i^T S_{ij} + S_{ji}^T F_j) + R_{ij}
                    matTemp5 = F_i'*Q_ij*F_j + (F_i'*S_ij+S_ji'*F_j) + R_ij;
                    lmiterm([-1, 2, n_j+1, 0],matTemp5);
                    
                    
                    % Filling W_ij^T = [W_ij1', W_ij3'; W_ij2', W_ij4']
                    % W_ij1' = - A_{ij}^T P_i -P_j B_j L_{ji}^T -L_ij^T + (-A_{ji}^T P_j + C_i^\T Q_{ij} C_j)^T;
                    lmiterm([-1, n_j, 1, P],-A_ij',1);
                    lmiterm([-1, n_j, 1, 0],matTemp1');
                    % lmiterm([-1, n_j, 1, -L_ji],-matTemp2',1);
                    % lmiterm([-1, 1, n_j, -L_ij],-1,1);
                    eval(['lmiterm([-1, 1, n_j, -L_',num2str(j),'i],-transpose(matTemp2),1);']);
                    eval(['lmiterm([-1, 1, n_j, -L_i',num2str(j),'],-1,1);']);
                    % W_ij3' = (-E_{ji}^\T P_j + S_{ji}^\T C_j + F_i^\T Q_{ji}^\T C_j)^T;
                    lmiterm([-1, n_j, 2, 0],matTemp4');
                    % W_ij2' = -E_{ij}^T P_i  + (C_i^\T S_{ij} + C_i^\T Q_{ij} F_j)^T;
                    lmiterm([-1, n_j+1, 1, P],-E_ij',1); 
                    lmiterm([-1, n_j+1, 1, 0],matTemp3'); 
                    % W_ij4' = (F_i^\T Q_{ij} F_j + (F_i^\T S_{ij} + S_{ji}^\T F_j) + R_{ij})^T
                    lmiterm([-1, n_j+1, 2, 0],matTemp5');
                    
                end

                % P term
                lmiterm([-2, 1, 1, P],1,1); % defines -P<0
                
                % Completing the LMI
                lmisys = getlmis;
                [tmin,sol] = feasp(lmisys); % Solve the LMI system
                isDissipative = tmin <= 0; % strictly feasible if this is satisfied
                P_i = dec2mat(lmisys,sol,P); % This needs to be stored
                L_i = dec2mat(lmisys,sol,L); % L_i = P_i B_i K_ii
                K_ii = pinv(B_i)*inv(P_i)*L_i;
                obj.controllerGains.localSFBDissLMI{iInd} = K_ii;
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    % L_ij = P_i B_i K_ij
                    % L_ij = dec2mat(lmisys,sol,L_ij); 
                    eval(['L_ij = dec2mat(lmisys,sol,L_i',num2str(j),');']);
                    K_ij = pinv(B_i)*inv(P_i)*L_ij;
                    K_ijVals{jInd} = K_ij;
                    L_ijVals{jInd} = L_ij;
                    obj.controllerGains.localSFBDissLMI{jInd} = K_ij;
                    
                    % L_ji = K_ji^T
                    % L_ji = dec2mat(lmisys,sol,L_ji); 
                    eval(['L_ji = dec2mat(lmisys,sol,L_',num2str(j),'i);']);
                    K_ji = transpose(L_ji);
                    K_jiVals{jInd} = K_ji; % these values will be loaded outside the function
                    L_jiVals{jInd} = L_ji;
                end
                % LMI End
                
                % Need to compute \tilede{W}_i and \tilde{W}_{ii} for storage
                
                % W_i term
                W_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % W_ij term
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};
                    C_j = subsystems(jInd).C{jInd};
                    E_ij = obj.E{jInd};
                    E_ji = subsystems(jInd).E{iInd};
                    F_j = subsystems(jInd).F{jInd};
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    Q_ij = obj.dataToBeDistributed.Q{jInd};
                    Q_ji = subsystems(jInd).dataToBeDistributed.Q{iInd};
                    S_ij = obj.dataToBeDistributed.S{jInd};
                    S_ji = subsystems(jInd).dataToBeDistributed.S{iInd};
                    R_ij = obj.dataToBeDistributed.R{jInd};
                    K_ij = K_ijVals{jInd}; % same as obj.localStabilizingSFBControllerGains{jInd};
                    K_ji = K_jiVals{jInd};
                    L_ij = L_ijVals{jInd};
                    L_ji = L_jiVals{jInd};
                    
                    % Filling W_ij = [W_ij1, W_ij2; W_ij3, W_ij4]
                    % W_ij1 = -A_{ji}^T P_j -L_{ji} B_j^T P_j - P_iA_{ij} - L_{ij} + C_i^T Q_{ij} C_j
                    W_ij1 = -P_i*A_ij -A_ji'*P_j - L_ji*B_j'*P_j -L_ij + C_i'*Q_ij*C_j;
                    % W_ij2 = -P_i E_{ij} + C_i^\T S_{ij} + C_i^\T Q_{ij} F_j;
                    W_ij2 = -P_i*E_ij + C_i'*S_ij + C_i'*Q_ij*F_j;
                    % W_ij3 = -E_{ji}^\T P_j + S_{ji}^\T C_j + F_i^\T Q_{ji}^\T C_j;
                    W_ij3 = -E_ji'*P_j + S_ji'*C_j + F_i'*Q_ji'*C_j;
                    % W_ij4 = F_i^\T Q_{ij} F_j + (F_i^\T S_{ij} + S_{ji}^\T F_j) + R_{ij};
                    W_ij4 = F_i'*Q_ij*F_j + (F_i'*S_ij+S_ji'*F_j) + R_ij;
                    
                    W_ij = [W_ij1, W_ij2; W_ij3, W_ij4];
                    W_i = [W_i, W_ij];
                end
                tildeW_i = W_i*M1_i;
                
                
                % W_ii term
                % W_ii1 = -A_{ii}^T P_i -K_{ii}^T B_i^T P_i - P_iA_{ii} - P_i B_i K_{ii} + C_i^T Q_{ii} C_i
                W_ii1 = -A_ii'*P_i - P_i*A_ii + C_i'*Q_ii*C_i - L_i - L_i'; % L_i = P_i B_i K_{ii}
                % W_ii2 = -P_i E_{ii} + C_i^T S_{ii} + C_i^T Q_{ii} F_i
                W_ii2 = -P_i*E_ii + C_i'*S_ii + C_i'*Q_ii*F_i;
                % W_ii3 = -E_{ii}^T P_i + S_{ii}^T C_i + F_i^T Q_{ii}^T C_i 
                W_ii3 = W_ii2';
                % W_ii4 = F_i^T Q_{ii} F_i + (F_i^T S_{ii} + S_{ii}^T F_i) + R_{ii}
                W_ii4 = F_i'*Q_ii*F_i + (F_i'*S_ii+S_ii'*F_i) + R_ii;
                
                W_ii = [W_ii1, W_ii2; W_ii3, W_ii4];
                
                tildeW_ii = W_ii - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                tildeW_i = [tildeW_i, tildeW_ii];

                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                if ~isDissipative
                    disp(['Not dissipative at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
                
            end
            
        end
            
            
        %% Local dissipative controller design using LMIs : an alternative approach using YALMIP
        function [isDissipative,K_ii,K_ijVals,K_jiVals] = designLocalDissipatingSFBControllers2(obj,previousSubsystems, subsystems)
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Dissipating at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            
            A_ii = obj.A{iInd};
            B_i = obj.B{iInd};
            C_i = obj.C{iInd};
            E_ii = obj.E{iInd};
            F_i = obj.F{iInd};
            Q_ii = obj.dataToBeDistributed.Q{iInd};
            S_ii = obj.dataToBeDistributed.Q{iInd};
            R_ii = obj.dataToBeDistributed.Q{iInd};
            
            dim_n = obj.dim_n;
            dim_p = obj.dim_p;
            dim_q = obj.dim_q;
            
            
            % For Dissipativity: W_ij = [W_ij1, W_ij2; W_ij3, W_ij4]
            % W_ij1 = -A_{ji}^T P_j -K_{ji}^T B_j^T P_j - P_iA_{ij} - P_i B_i K_{ij} + C_i^T Q_{ij} C_j
            % W_ij2 = -P_i E_{ij} + C_i^T S_{ij} + C_i^T Q_{ij} F_j
            % W_ij3 = -E_{ji}^T P_j + S_{ji}^T C_j + F_i^T Q_{ji}^T C_j  
            % W_ij4 = F_i^T Q_{ij} F_j + (F_i^T S_{ij} + S_{ji}^T F_j) + R_{ij}
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii > 0 and P_i > 0 with K_ii free
                % W_ii = [W_ii1, W_ii2; W_ii3, W_ii4]
                % W_ii1 = -A_{ii}^T P_i -K_{ii}^T B_i^T P_i - P_iA_{ii} - P_i B_i K_{ii} + C_i^T Q_{ii} C_i
                % W_ii2 = -P_i E_{ii} + C_i^T S_{ii} + C_i^T Q_{ii} F_i
                % W_ii3 = -E_{ii}^T P_i + S_{ii}^T C_i + F_i^T Q_{ii}^T C_i  
                % W_ii4 = F_i^T Q_{ii} F_i + (F_i^T S_{ii} + S_{ii}^T F_i) + R_{ii}
            
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                % P = lmivar(1,[dim_n, 1]);   % P is a variable, 1: square symmetric, size(A,1) gives the size and 1 gives that Q is a full matrix 
                P = sdpvar(dim_n,dim_n)
                
                
                
                
                
                
                
                
                
                
                
                L = lmivar(2,[dim_n,dim_n])      % L_i = P_i B_i K_ii is a variable, rectangular matrix
                
                
                % W_ii1 = -A_{ii}^T P_i -K_{ii}^T B_i^T P_i - P_iA_{ii} - P_i B_i K_{ii} + C_i^T Q_{ii} C_i
                % W_ii1 = -A_{ii}^T P_i -L_i^T - P_iA_{ii} - L_i + C_i^T Q_{ii} C_i
                lmiterm([-1, 1, 1, P],-1,A_ii,'s');
                lmiterm([-1, 1, 1, L],-1,1,'s');
                matTemp1 = C_i'*Q_ii*C_i;
                lmiterm([-1, 1, 1, 0],matTemp1);
                % W_ii2 = -P_i E_{ii} + C_i^T S_{ii} + C_i^T Q_{ii} F_i
                lmiterm([-1, 1, 2, P],-1,E_ii); 
                matTemp2 = C_i'*S_ii + C_i'*Q_ii*F_i;
                lmiterm([-1, 1, 2, 0],matTemp2); 
                % W_ii3 = -E_{ii}^T P_i + S_{ii}^T C_i + F_i^T Q_{ii}^T C_i  
                lmiterm([-1, 2, 1, P],-E_ii',1);
                lmiterm([-1, 2, 1, 0],matTemp2'); 
                % W_ii4 = F_i^T Q_{ii} F_i + (F_i^T S_{ii} + S_{ii}^T F_i) + R_{ii}
                matTemp4 = F_i'*Q_ii*F_i + (F_i'*S_ii+S_ii'*F_i) + R_ii;
                lmiterm([-1, 2, 2, 0],matTemp4);
                
                % P>0
                lmiterm([-2, 1, 1, P],1,1); % defines -P<0
                lmisys = getlmis;
                [tmin,sol] = feasp(lmisys); % Solve the LMI system
                isDissipative = tmin <= 0; % strictly feasible if this is satisfied
                
                P_i = dec2mat(lmisys,sol,P); % This needs to be stored
                L_i = dec2mat(lmisys,sol,L); % This needs to be stored
                K_ii = pinv(B_i)*inv(P_i)*L_i;
                K_jiVals = [];
                K_ijVals = [];
                tildeW_ii1 = -A_ii'*P_i - P_i*A_ii + L_i + L_i' + matTemp1;
                tildeW_ii2 = -P_i*E_ii + matTemp2;
                tildeW_ii3 = tildeW_ii2';
                tildeW_ii4 = matTemp4;
                tildeW_i = [tildeW_ii1, tildeW_ii2; tildeW_ii3, tildeW_ii4]; % Note that here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                % LMI End
                
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                obj.controllerGains.localSFBDissLMI{iInd} = K_ii;
                if ~isDissipative
                    disp(['Not dissipative at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
                   
            else
                % This subsystem has to talk with all the previosSubsystems
                % tildeW_ii > 0 iff [W_ii, W_i; W_i', M_i] > 0 is required where 
                % M_i = inv((scriptD_i*scriptA_i)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i)^{-1}')
                
                % Introduce L_i = P_i B_i K_ii;  L_ij = P_i B_i K_ij; L_ji = K_ji' as variables apart from variables P_i.
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                P = lmivar(1,[dim_n, 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    % L_ij = lmivar(2,[dim_n, dim_n])     % For L_ij = P_i B_i K_ij
                    % L_ji = lmivar(2,[dim_n, dim_p])     % For L_ji = K_ji^T
                    eval(['L_i',num2str(j),'= lmivar(2,[dim_n, dim_n]);'])
                    eval(['L_',num2str(j),'i= lmivar(2,[dim_n, dim_p]);'])
                end
                L = lmivar(2,[dim_n,dim_n]); % For L_i = P_i B_i K_ii 
                                
                % W_ii term: W_ii = [W_ii1, W_ii2; W_ii3, W_ii4]
                % W_ii1 = -A_{ii}^T P_i -K_{ii}^T B_i^T P_i - P_iA_{ii} - P_i B_i K_{ii} + C_i^T Q_{ii} C_i
                % W_ii1 = -A_{ii}^T P_i -L_i^T - P_iA_{ii} - L_i + C_i^T Q_{ii} C_i
                lmiterm([-1, 1, 1, P],-1,A_ii,'s');
                lmiterm([-1, 1, 1, L],-1,1,'s');
                matTemp1 = C_i'*Q_ii*C_i;
                lmiterm([-1, 1, 1, 0],matTemp1);
                % W_ii2 = -P_i E_{ii} + C_i^T S_{ii} + C_i^T Q_{ii} F_i
                lmiterm([-1, 1, 2, P],-1,E_ii); 
                matTemp2 = C_i'*S_ii + C_i'*Q_ii*F_i;
                lmiterm([-1, 1, 2, 0],matTemp2); 
                % W_ii3 = -E_{ii}^T P_i + S_{ii}^T C_i + F_i^T Q_{ii}^T C_i  
                lmiterm([-1, 2, 1, P],-E_ii',1);
                lmiterm([-1, 2, 1, 0],matTemp2'); 
                % W_ii4 = F_i^T Q_{ii} F_i + (F_i^T S_{ii} + S_{ii}^T F_i) + R_{ii}
                matTemp4 = F_i'*Q_ii*F_i + (F_i'*S_ii+S_ii'*F_i) + R_ii;
                lmiterm([-1, 2, 2, 0],matTemp4);
                
                
                % M_i term = inv((scriptD_i*scriptA_i)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i)^{-1}')
                blockSize = dim_n+dim_q;
                scriptA_i = [];
                scriptD_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % Getting stored info from jInd to create \mathcal{A}_i and \mathcal{D}_i matrices (their j-th columns)
                    tildeW_j = subsystems(jInd).dataToBeDistributed.tildeW;
                                        
                    Z = zeros(blockSize*(i-1-j),blockSize); % (i-1)-j blocks of blockSizeXblockSize zero matrices
                    z = zeros(blockSize*(j-1),blockSize);
                    if j==1
                        tildeW_jj = tildeW_j;                    
                        scriptA_i = [tildeW_jj; Z];         % The first column of \mathcal{A}_i.
                        scriptD_i = [inv(tildeW_jj); Z];    % The first column of \mathcal{D}_i.
                    else
                        tildeW_jj = tildeW_j(:,blockSize*(j-1)+1:blockSize*j);   % last blockSizeXblockSize block in the row block vector
                        tildeW_j  = tildeW_j(:,1:blockSize*(j-1));                % first (j-1) blocks in the row block vector
                        scriptA_i = [scriptA_i, [tildeW_j'; tildeW_jj ; Z]];    % The j-th column of \mathcal{A}_i.
                        scriptD_i = [scriptD_i, [z; inv(tildeW_jj); Z]];         % The j-th column of \mathcal{D}_i.
                    end                    
                end
                M1_i = inv(scriptD_i*scriptA_i);
                % M_i = inv(M1_i*scriptD_i*M1_i'); % THis fills 2*(i-1)x2*(i-1) blocks in the LMI
                M_i = scriptA_i'*scriptD_i'*scriptA_i;
                if issymmetric(scriptD_i) & issymmetric(scriptA_i) & ~issymmetric(M_i)
                    tf = norm(M_i-M_i.',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M_i = 0.5*(M_i + M_i');
                end
                
                % Filling M_i
                for x=1:1:length(previousSubsystems)
                    m_x = (x-1)*blockSize; % exact index
                    n_x = 2+2*(x-1)+1;           % block index
                    for y=1:1:length(previousSubsystems)
                        m_y = (y-1)*blockSize;  % exact index
                        n_y = 2+2*(y-1)+1;           % block index
                        M_i_xy = M_i(m_x+1:m_x+blockSize, m_y+1:m_y+blockSize);
                        M_i_xy1 = M_i_xy(1:obj.dim_n,1:obj.dim_n);
                        M_i_xy2 = M_i_xy(1:obj.dim_n,obj.dim_n+1:blockSize);
                        M_i_xy3 = M_i_xy(obj.dim_n+1:blockSize,1:obj.dim_n);
                        M_i_xy4 = M_i_xy(obj.dim_n+1:blockSize,obj.dim_n+1:blockSize);
                        lmiterm([-1, n_x, n_y, 0],M_i_xy1);
                        lmiterm([-1, n_x, n_y+1, 0],M_i_xy2);
                        lmiterm([-1, n_x+1, n_y, 0],M_i_xy3);
                        lmiterm([-1, n_x+1, n_y+1, 0],M_i_xy4);
                    end
                end
               
                % W_i term
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % W_ij term
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};
                    B_j = subsystems(jInd).B{jInd};
                    C_j = subsystems(jInd).C{jInd};
                    E_ij = obj.E{jInd};
                    E_ji = subsystems(jInd).E{iInd};
                    F_j = subsystems(jInd).F{jInd};
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    Q_ij = obj.dataToBeDistributed.Q{jInd};
                    Q_ji = subsystems(jInd).dataToBeDistributed.Q{iInd};
                    S_ij = obj.dataToBeDistributed.S{jInd};
                    S_ji = subsystems(jInd).dataToBeDistributed.S{iInd};
                    R_ij = obj.dataToBeDistributed.R{jInd};
                    n_j = 2+2*(j-1)+1;
                    
                    
                    % Filling W_ij = [W_ij1, W_ij2; W_ij3, W_ij4]
                    % W_ij1 = -A_{ji}^T P_j -K_{ji}^T B_j^T P_j - P_iA_{ij} - P_i B_i K_{ij} + C_i^T Q_{ij} C_j
                    % W_ij1 = -A_{ji}^T P_j -L_{ji} B_j^T P_j - P_iA_{ij} - L_{ij} + C_i^T Q_{ij} C_j
                    lmiterm([-1, 1, n_j, P],-1,A_ij);
                    matTemp1 = -A_ji'*P_j + C_i'*Q_ij*C_j;
                    lmiterm([-1, 1, n_j, 0],matTemp1);
                    matTemp2 = B_j'*P_j;
                    % lmiterm([-1, 1, n_j, L_ji],-1,matTemp2);
                    % lmiterm([-1, 1, n_j, L_ij],-1,1);
                    eval(['lmiterm([-1, 1, n_j, L_',num2str(j),'i],-1,matTemp2);']);
                    eval(['lmiterm([-1, 1, n_j, L_i',num2str(j),'],-1,1);']);
                    % W_ij2 = -P_i E_{ij} + C_i^T S_{ij} + C_i^T Q_{ij} F_j
                    lmiterm([-1, 1, n_j+1, P],-1,E_ij); 
                    matTemp3 = C_i'*S_ij + C_i'*Q_ij*F_j;
                    lmiterm([-1, 1, n_j+1, 0],matTemp3);
                    % W_ij3 = -E_{ji}^T P_j + S_{ji}^T C_j + F_i^T Q_{ji}^T C_j  
                    matTemp4 = -E_ji'*P_j + S_ji'*C_j + F_i'*Q_ji'*C_j;
                    lmiterm([-1, 2, n_j, 0],matTemp4);                    
                    % W_ij4 = F_i^T Q_{ij} F_j + (F_i^T S_{ij} + S_{ji}^T F_j) + R_{ij}
                    matTemp5 = F_i'*Q_ij*F_j + (F_i'*S_ij+S_ji'*F_j) + R_ij;
                    lmiterm([-1, 2, n_j+1, 0],matTemp5);
                    
                    
                    % Filling W_ij^T = [W_ij1', W_ij3'; W_ij2', W_ij4']
                    % W_ij1' = - A_{ij}^T P_i -P_j B_j L_{ji}^T -L_ij^T + (-A_{ji}^T P_j + C_i^\T Q_{ij} C_j)^T;
                    lmiterm([-1, n_j, 1, P],-A_ij',1);
                    lmiterm([-1, n_j, 1, 0],matTemp1');
                    % lmiterm([-1, n_j, 1, -L_ji],-matTemp2',1);
                    % lmiterm([-1, 1, n_j, -L_ij],-1,1);
                    eval(['lmiterm([-1, 1, n_j, -L_',num2str(j),'i],-transpose(matTemp2),1);']);
                    eval(['lmiterm([-1, 1, n_j, -L_i',num2str(j),'],-1,1);']);
                    % W_ij3' = (-E_{ji}^\T P_j + S_{ji}^\T C_j + F_i^\T Q_{ji}^\T C_j)^T;
                    lmiterm([-1, n_j, 2, 0],matTemp4');
                    % W_ij2' = -E_{ij}^T P_i  + (C_i^\T S_{ij} + C_i^\T Q_{ij} F_j)^T;
                    lmiterm([-1, n_j+1, 1, P],-E_ij',1); 
                    lmiterm([-1, n_j+1, 1, 0],matTemp3'); 
                    % W_ij4' = (F_i^\T Q_{ij} F_j + (F_i^\T S_{ij} + S_{ji}^\T F_j) + R_{ij})^T
                    lmiterm([-1, n_j+1, 2, 0],matTemp5');
                    
                end

                % P term
                lmiterm([-2, 1, 1, P],1,1); % defines -P<0
                
                % Completing the LMI
                lmisys = getlmis;
                [tmin,sol] = feasp(lmisys); % Solve the LMI system
                isDissipative = tmin <= 0; % strictly feasible if this is satisfied
                P_i = dec2mat(lmisys,sol,P); % This needs to be stored
                L_i = dec2mat(lmisys,sol,L); % L_i = P_i B_i K_ii
                K_ii = pinv(B_i)*inv(P_i)*L_i;
                obj.controllerGains.localSFBDissLMI{iInd} = K_ii;
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    % L_ij = P_i B_i K_ij
                    % L_ij = dec2mat(lmisys,sol,L_ij); 
                    eval(['L_ij = dec2mat(lmisys,sol,L_i',num2str(j),');']);
                    K_ij = pinv(B_i)*inv(P_i)*L_ij;
                    K_ijVals{jInd} = K_ij;
                    L_ijVals{jInd} = L_ij;
                    obj.controllerGains.localSFBDissLMI{jInd} = K_ij;
                    
                    % L_ji = K_ji^T
                    % L_ji = dec2mat(lmisys,sol,L_ji); 
                    eval(['L_ji = dec2mat(lmisys,sol,L_',num2str(j),'i);']);
                    K_ji = transpose(L_ji);
                    K_jiVals{jInd} = K_ji; % these values will be loaded outside the function
                    L_jiVals{jInd} = L_ji;
                end
                % LMI End
                
                % Need to compute \tilede{W}_i and \tilde{W}_{ii} for storage
                
                % W_i term
                W_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % W_ij term
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};
                    C_j = subsystems(jInd).C{jInd};
                    E_ij = obj.E{jInd};
                    E_ji = subsystems(jInd).E{iInd};
                    F_j = subsystems(jInd).F{jInd};
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    Q_ij = obj.dataToBeDistributed.Q{jInd};
                    Q_ji = subsystems(jInd).dataToBeDistributed.Q{iInd};
                    S_ij = obj.dataToBeDistributed.S{jInd};
                    S_ji = subsystems(jInd).dataToBeDistributed.S{iInd};
                    R_ij = obj.dataToBeDistributed.R{jInd};
                    K_ij = K_ijVals{jInd}; % same as obj.localStabilizingSFBControllerGains{jInd};
                    K_ji = K_jiVals{jInd};
                    L_ij = L_ijVals{jInd};
                    L_ji = L_jiVals{jInd};
                    
                    % Filling W_ij = [W_ij1, W_ij2; W_ij3, W_ij4]
                    % W_ij1 = -A_{ji}^T P_j -L_{ji} B_j^T P_j - P_iA_{ij} - L_{ij} + C_i^T Q_{ij} C_j
                    W_ij1 = -P_i*A_ij -A_ji'*P_j - L_ji*B_j'*P_j -L_ij + C_i'*Q_ij*C_j;
                    % W_ij2 = -P_i E_{ij} + C_i^\T S_{ij} + C_i^\T Q_{ij} F_j;
                    W_ij2 = -P_i*E_ij + C_i'*S_ij + C_i'*Q_ij*F_j;
                    % W_ij3 = -E_{ji}^\T P_j + S_{ji}^\T C_j + F_i^\T Q_{ji}^\T C_j;
                    W_ij3 = -E_ji'*P_j + S_ji'*C_j + F_i'*Q_ji'*C_j;
                    % W_ij4 = F_i^\T Q_{ij} F_j + (F_i^\T S_{ij} + S_{ji}^\T F_j) + R_{ij};
                    W_ij4 = F_i'*Q_ij*F_j + (F_i'*S_ij+S_ji'*F_j) + R_ij;
                    
                    W_ij = [W_ij1, W_ij2; W_ij3, W_ij4];
                    W_i = [W_i, W_ij];
                end
                tildeW_i = W_i*M1_i;
                
                
                % W_ii term
                % W_ii1 = -A_{ii}^T P_i -K_{ii}^T B_i^T P_i - P_iA_{ii} - P_i B_i K_{ii} + C_i^T Q_{ii} C_i
                W_ii1 = -A_ii'*P_i - P_i*A_ii + C_i'*Q_ii*C_i - L_i - L_i'; % L_i = P_i B_i K_{ii}
                % W_ii2 = -P_i E_{ii} + C_i^T S_{ii} + C_i^T Q_{ii} F_i
                W_ii2 = -P_i*E_ii + C_i'*S_ii + C_i'*Q_ii*F_i;
                % W_ii3 = -E_{ii}^T P_i + S_{ii}^T C_i + F_i^T Q_{ii}^T C_i 
                W_ii3 = W_ii2';
                % W_ii4 = F_i^T Q_{ii} F_i + (F_i^T S_{ii} + S_{ii}^T F_i) + R_{ii}
                W_ii4 = F_i'*Q_ii*F_i + (F_i'*S_ii+S_ii'*F_i) + R_ii;
                
                W_ii = [W_ii1, W_ii2; W_ii3, W_ii4];
                
                tildeW_ii = W_ii - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                tildeW_i = [tildeW_i, tildeW_ii];

                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                if ~isDissipative
                    disp(['Not dissipative at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
                
            end
            
        end
            
          
        
        
        
        
    end
end

