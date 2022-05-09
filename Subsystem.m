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
        controllerGains
        
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
                        sys1 = rss(n,m,p);
                        while ~all(abs(real(eig(sys1.A)))>0.1)
                            sys1 = rss(n,m,p);
                        end
                        sys2 = rss(n,m,q);
                        obj.A{j} = sys1.A; % A_ii
                        obj.B{j} = sys1.B;
                        obj.C{j} = sys1.C;
                        obj.D{j} = sys1.D;
                        obj.E{j} = 0.5*sys2.B;
                        obj.F{j} = 0.1*sys2.D;
                    else
                        [n_j, p_j, q_j, m_j] = subsystems(j).getDimensions();
                        sys1 = rss(n,m_j,p_j);  % assume n=n_j
                        while ~all(abs(real(eig(sys1.A)))>0.1)
                            sys1 = rss(n,m_j,p_j);  % assume n=n_j
                        end
                        sys2 = rss(n,m_j,q_j); % assume n=n_j
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
                        [n_j, p_j, q_j, m_j] = subsystems(j).getDimensions();
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
            obj.controllerGains.localSFBLQR = K;
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
            text(obj.location(1)-0.035,obj.location(2)-0.03,num2str(obj.index),'Color','k','FontSize',10);
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
        
        
        %% An old local stability check approach
        function isFeasible = checkStability(obj,previousSubsystems,subsystems)
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Checking stability at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            % For stability: W_ij = -A_ji^T P_j -P_i A_ij
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii
                % W_ii = -A_ii^T P_i -P_i A_ii
                % we need M = W_ii>0 and P_i > 0
                A_ii = obj.A{iInd};
                
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                P = lmivar(1,[size(A_ii,1), 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 
                
                lmiterm([-1, 1, 1, P],-1,A_ii,'s'); % defines -(-PA-A'P)<0; here "s" flag is used to get the symmetric part of PA.
                lmiterm([-2, 1, 1, P],1,1); % defines -P<0
                lmisys = getlmis;
                
                [tmin,Psol] = feasp(lmisys); % Solve the LMI system
                isFeasible = tmin <= 0; % strictly feasible if this is satisfied
                P_i = dec2mat(lmisys,Psol,P); % This needs to be stored
                tildeW_i = -A_ii'*P_i - P_i*A_ii; % Note that here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                % LMI End
                
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                if ~isFeasible
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
                   
            else
                % This subsystem has to talk with all the previosSubsystems
                scriptA_i = [];
                scriptD_i = [];
                W_i_1 = []; % W_i^T = W_i_1 + P_i X W_i_2 (both W_i_1 and W_i_2 are row block vectors, X is kronecker block product)
                W_i_2 = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % Getting stored info from jInd to create \mathcal{A}_i and \mathcal{D}_i matrices (their j-th columns)
                    tildeW_j = subsystems(jInd).dataToBeDistributed.tildeW;
                    Z = zeros(obj.dim_n*(i-1-j),obj.dim_n); % (i-1)-j blocks of nxn zero matrices
                    z = zeros(obj.dim_n*(j-1),obj.dim_n);
                    if j==1
                        tildeW_jj = tildeW_j;                    
                        scriptA_i = [tildeW_jj; Z];         % The first column of \mathcal{A}_i.
                        scriptD_i = [inv(tildeW_jj); Z];    % The first column of \mathcal{D}_i.
                    else
                        tildeW_jj = tildeW_j(:,obj.dim_n*(j-1)+1:obj.dim_n*j);   % last nxn block in the row block vector
                        tildeW_j  = tildeW_j(:,1:obj.dim_n*(j-1));                % first (j-1) blocks in the row block vector
                        scriptA_i = [scriptA_i, [tildeW_j'; tildeW_jj ; Z]];    % The j-th column of \mathcal{A}_i.
                        scriptD_i = [scriptD_i, [z; inv(tildeW_jj); Z]];         % The j-th column of \mathcal{D}_i.
                    end

                    % Getting required information from jInd to compute W_ij (to eventually create W_i) 
                    % Note that: W_ij = -A_ji^T P_j - P_i A_ij
                    A_ji = subsystems(jInd).A{iInd};
                    A_ij = obj.A{jInd};
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    W_ij_1 = -A_ji'*P_j;
                    W_i_1 = [W_i_1, W_ij_1];
                    W_ij_2 = -A_ij; % note that we cannot pre-multiply this by P_i as P_i is a variable here.
                    W_i_2 = [W_i_2, W_ij_2];
                end

                % So far, we have W_i (its components) and \mathcal{A}_i and \mathcal{D}_i matrices. 
                % W_ii can be computed easily using: W_ii = -A_ii^T P_i - P_i A_ii (will be a function of P_i)
                A_ii = obj.A{iInd};
                % Since W_i = W_i_1 + P X W_i_2 and \tildeW_i = W_i*(\mathcal{D}_i \mathcal{A}_i)^{-1}
                % and \tilde{W}_ii = W_ii - \tilde{W}_i \mathcal{D}_i \tilde{W}_i^T we get:
                M_i = inv(scriptD_i*scriptA_i);
                Mat1 = -W_i_1*M_i*scriptD_i*M_i'*W_i_1';
                Mat2 = -(A_ii' + W_i_1*M_i*scriptD_i*M_i'*W_i_2');
                Mat3 = -(A_ii + W_i_2*M_i*scriptD_i*M_i'*W_i_1');
                Mat4 = W_i_2*M_i;
                Mat5 = -scriptD_i;

%                 disp("Interested matrices: ")
%                 scriptA_i
%                 scriptD_i
%                 M_i
%                 W_i_1
%                 W_i_2
                
                
                % Interested LMI: P>0 such that Mat1 +  Mat2*P + P*Mat3 + P*Mat4*Mat5*Mat4'*P > 0 (need to use the schur complement)
                % Otherwise: P>0 such that [Mat1 +  Mat2*P + P*Mat3,  P*Mat4;  Mat4'*P,   -inv(Mat5) ]
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                P = lmivar(1,[size(A_ii,1), 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 

                lmiterm([-1, 1, 1, P],Mat2,1,'s'); % defines Mat2*P + P*Mat3
                lmiterm([-1, 1, 1, 0],Mat1); % defines Mat1
                lmiterm([-1, 1, 2, P],1,Mat4);
                lmiterm([-1, 2, 1, P],Mat4',1);
                lmiterm([-1, 2, 2, 0],-inv(Mat5));
                lmiterm([-2, 1, 1, P],1,1); % defines -P<0
                lmisys = getlmis;

                [tmin,Psol] = feasp(lmisys); % Solve the LMI system
                isFeasible = tmin <= 0; % strictly feasible if this is satisfied
                P_i = dec2mat(lmisys,Psol,P); % This needs to be stored
                % LMI End

                % Need to compute \tilede{W}_i and \tilde{W}_{ii} for storage
                % P_i X W_i_2
                W_i_3 = [];
                for j = 1:1:length(previousSubsystems)
                    W_i_3 = [W_i_3, P_i*W_i_2(:,obj.dim_n*(j-1)+1:obj.dim_n*j)];
                end
                W_i = W_i_1 + W_i_3;
                tildeW_i = W_i*M_i;
                W_ii = -A_ii'*P_i - P_i*A_ii;
                tildeW_ii = W_ii - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                tildeW_i = [tildeW_i, tildeW_ii];

                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                if ~isFeasible
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            end
          
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
        
        %% Local stability check: LMI approach
        % The alternative approach to check the stability (via schur's compliment)
        function isFeasible = checkStability2(obj,previousSubsystems,subsystems)
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Checking stability at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            A_ii = obj.A{iInd};
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii
                % W_ii = -A_ii^T P_i -P_i A_ii
                % we need M = W_ii>0 and P_i > 0
                
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                P = lmivar(1,[size(A_ii,1), 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 
                
                lmiterm([-1, 1, 1, P],-1,A_ii,'s'); % defines -(-PA-A'P)<0; here "s" flag is used to get the symmetric part of PA.
                lmiterm([-2, 1, 1, P],1,1); % defines -P<0
                lmisys = getlmis;
                
                [tmin,Psol] = feasp(lmisys); % Solve the LMI system
                isFeasible = tmin <= 0; % strictly feasible if this is satisfied
                P_i = dec2mat(lmisys,Psol,P); % This needs to be stored
                disp(['Data saved ',num2str(iInd),' after ',num2str(previousSubsystems),'.'])
                tildeW_i = -A_ii'*P_i - P_i*A_ii % Note that here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                % LMI End
                
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                if ~isFeasible
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            
            else
                % This subsystem has to talk with all the previosSubsystems
                % tildeW_ii > 0 iff [W_ii, W_i; W_i', M_i] > 0 is required where 
                % M_i = inv((scriptD_i*scriptA_i)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i)^{-1}')
                
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                P = lmivar(1,[size(A_ii,1), 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 
                
                % W_ii term
                lmiterm([-1, 1, 1, P],-1,A_ii,'s'); % defines -(-PA-A'P)<0; here "s" flag is used to get the symmetric part of PA.
                
                % M_i term
                blockSize = obj.dim_n; 
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
                disp(['Data at ',num2str(iInd),' after ',num2str(previousSubsystems),'.'])
                scriptA_i;
                scriptD_i;
                
                M1_i = inv(scriptD_i*scriptA_i);
                % M_i = inv(M1_i*scriptD_i*M1_i') % THis fills (i-1)x(i-1) blocks in the LMI
                M_i = scriptA_i'*scriptD_i*scriptA_i;
                if issymmetric(scriptD_i) & issymmetric(scriptA_i) & ~issymmetric(M_i)
                    tf = norm(M_i-M_i.',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M_i = 0.5*(M_i + M_i');
                end
                % save('test','M_i','scriptD_i','scriptA_i')
                
                % Filling M_i
                for x=1:1:length(previousSubsystems)
                    m_x = (x-1)*blockSize; % exact index
                    n_x = 1+(x-1)+1;           % block index
                    for y=1:1:length(previousSubsystems)
                        m_y = (y-1)*blockSize;  % exact index
                        n_y = 1+(y-1)+1;           % block index
                        M_i_xy = M_i(m_x+1:m_x+blockSize, m_y+1:m_y+blockSize);
                        lmiterm([-1, n_x, n_y, 0],M_i_xy);
                    end
                end
                
                
                % W_i term
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % W_ij term
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};                    
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    n_j = 1+(j-1)+1;
                    
                    % Filling W_ij = -A_{ji}^T P_j - P_iA_{ij}
                    lmiterm([-1, 1, n_j, P],-1,A_ij);
                    matTemp1 = -A_ji'*P_j;
                    lmiterm([-1, 1, n_j, 0],matTemp1);
                    
                    % Filling W_ij^T = (-A_{ji}^T P_j - P_iA_{ij})^T
                    lmiterm([-1, n_j, 1, P],-A_ij',1);
                    lmiterm([-1, n_j, 1, 0],matTemp1');
                end
                
                % P term
                lmiterm([-2, 1, 1, P],1,1); % defines -P<0
                
                % Completing the LMI
                lmisys = getlmis;
%                 options = [0,0,100,0,0]; %omit,maxiter,rad,...
                [tmin,Psol] = feasp(lmisys); % Solve the LMI system
                isFeasible = tmin <= 0; % strictly feasible if this is satisfied
                P_i = dec2mat(lmisys,Psol,P); % This needs to be stored
                % LMI End
                
                % Need to compute \tilede{W}_i and \tilde{W}_{ii} for storage
                
                % W_i term
                W_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % W_ij term
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    
                    % Filling W_ij = -A_{ji}^T P_j - P_iA_{ij}
                    W_ij = -P_i*A_ij - A_ji'*P_j;
                    W_i = [W_i, W_ij];
                end
                tildeW_i = W_i*M1_i;
                
                % W_ii term: W_ii = -A_{ii}^T P_i - P_iA_{ii}
                W_ii = -A_ii'*P_i - P_i*A_ii;
                
                tildeW_ii = W_ii - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                disp(['Data saved ',num2str(iInd),' after ',num2str(previousSubsystems),'.'])
                tildeW_i = [tildeW_i, tildeW_ii]

                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                if ~isFeasible
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            end
            
            
            
        end
        
        
        %% Local QSR-Dissipativity check: LMI approach
        function isFeasible = checkQSRDissipativity(obj,previousSubsystems,subsystems)
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Checking QSR-Dissipativity at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            
            A_ii = obj.A{iInd};
            C_i = obj.C{iInd};
            E_ii = obj.E{iInd};
            F_i = obj.F{iInd};
            Q_ii = obj.dataToBeDistributed.Q{iInd};
            S_ii = obj.dataToBeDistributed.Q{iInd};
            R_ii = obj.dataToBeDistributed.Q{iInd};
            
            % For Dissipativity: W_ij = [W_ij1, W_ij2; W_ij3, W_ij4]
            % W_ij1 = -A_{ji}^T P_j - P_iA_{ij} + C_i^\T Q_{ij} C_j;
            % W_ij2 = -P_i E_{ij} + C_i^\T S_{ij} + C_i^\T Q_{ij} F_j;
            % W_ij3 = -E_{ji}^\T P_j + S_{ji}^\T C_j + F_i^\T Q_{ji}^\T C_j;
            % W_ij4 = F_i^\T Q_{ij} F_j + (F_i^\T S_{ij} + S_{ji}^\T F_j) + R_{ij};
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii > 0
                % W_ii = [W_ii1, W_ii2; W_ii3, W_ii4]
                % W_ii1 = -A_{ii}^T P_i - P_iA_{ii} + C_i^\T Q_{ii} C_i;
                % W_ii2 = -P_i E_{ii} + C_i^\T S_{ii} + C_i^\T Q_{ii} F_i;
                % W_ii3 = -E_{ii}^\T P_i + S_{ii}^\T C_i + F_i^\T Q_{ii}^\T C_i;
                % W_ii4 = F_i^\T Q_{ii} F_i + (F_i^\T S_{ii} + S_{ii}^\T F_i) + R_{ii};
                
                % we need M = W_ii>0 and P_i > 0
                
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                P = lmivar(1,[size(A_ii,1), 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 
                
                % W_ii1 = -A_{ii}^T P_i - P_iA_{ii} + C_i^\T Q_{ii} C_i
                lmiterm([-1, 1, 1, P],-1,A_ii,'s'); 
                matTemp1 = C_i'*Q_ii*C_i;
                lmiterm([-1, 1, 1, 0],matTemp1);
                % W_ii2 = -P_i E_{ii} + C_i^\T S_{ii} + C_i^\T Q_{ii} F_i;
                lmiterm([-1, 1, 2, P],-1,E_ii); 
                matTemp2 = C_i'*S_ii + C_i'*Q_ii*F_i;
                lmiterm([-1, 1, 2, 0],matTemp2); 
                % W_ii3 = -E_{ii}^\T P_i + S_{ii}^\T C_i + F_i^\T Q_{ii}^\T C_i;
                lmiterm([-1, 2, 1, P],-E_ii',1);
                lmiterm([-1, 2, 1, 0],matTemp2'); 
                % W_ii4 = F_i^\T Q_{ii} F_i + (F_i^\T S_{ii} + S_{ii}^\T F_i) + R_{ii};
                matTemp4 = F_i'*Q_ii*F_i + (F_i'*S_ii+S_ii'*F_i) + R_ii;
                lmiterm([-1, 2, 2, 0],matTemp4);
                
                % P>0
                lmiterm([-2, 1, 1, P],1,1); % defines -P<0
                lmisys = getlmis;
                [tmin,Psol] = feasp(lmisys); % Solve the LMI system
                isFeasible = tmin <= 0; % strictly feasible if this is satisfied
                
                P_i = dec2mat(lmisys,Psol,P); % This needs to be stored
                tildeW_ii1 = -A_ii'*P_i - P_i*A_ii + matTemp1;
                tildeW_ii2 = -P_i*E_ii + matTemp2;
                tildeW_ii3 = tildeW_ii2';
                tildeW_ii4 = matTemp4;
                tildeW_i = [tildeW_ii1, tildeW_ii2; tildeW_ii3, tildeW_ii4]; % Note that here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                % LMI End
                
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                if ~isFeasible
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
                   
            else
                disp("Here!!!")
                % This subsystem has to talk with all the previosSubsystems
                % tildeW_ii > 0 iff [W_ii, W_i; W_i', M_i] > 0 is required where 
                % M_i = inv((scriptD_i*scriptA_i)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i)^{-1}')
                
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                P = lmivar(1,[size(A_ii,1), 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 
                                
                % W_ii term
                % W_ii1 = -A_{ii}^T P_i - P_iA_{ii} + C_i^\T Q_{ii} C_i
                lmiterm([-1, 1, 1, P],-1,A_ii,'s'); 
                matTemp1 = C_i'*Q_ii*C_i;
                lmiterm([-1, 1, 1, 0],matTemp1);
                % W_ii2 = -P_i E_{ii} + C_i^\T S_{ii} + C_i^\T Q_{ii} F_i;
                lmiterm([-1, 1, 2, P],-1,E_ii); 
                matTemp2 = C_i'*S_ii + C_i'*Q_ii*F_i;
                lmiterm([-1, 1, 2, 0],matTemp2); 
                % W_ii3 = -E_{ii}^\T P_i + S_{ii}^\T C_i + F_i^\T Q_{ii}^\T C_i;
                lmiterm([-1, 2, 1, P],-E_ii',1);
                lmiterm([-1, 2, 1, 0],matTemp2'); 
                % W_ii4 = F_i^\T Q_{ii} F_i + (F_i^\T S_{ii} + S_{ii}^\T F_i) + R_{ii};
                matTemp4 = F_i'*Q_ii*F_i + (F_i'*S_ii+S_ii'*F_i) + R_ii;
                lmiterm([-1, 2, 2, 0],matTemp4);
                
                
                % M_i term = inv((scriptD_i*scriptA_i)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i)^{-1}')
                blockSize = obj.dim_n+obj.dim_q; 
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
%                 M_i = inv(M1_i*scriptD_i*M1_i'); % THis fills 2*(i-1)x2*(i-1) blocks in the LMI
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
                    % W_ij1 = -A_{ji}^T P_j - P_iA_{ij} + C_i^\T Q_{ij} C_j;
                    lmiterm([-1, 1, n_j, P],-1,A_ij);
                    matTemp1 = -A_ji'*P_j + C_i'*Q_ij*C_j;
                    lmiterm([-1, 1, n_j, 0],matTemp1);
                    % W_ij2 = -P_i E_{ij} + C_i^\T S_{ij} + C_i^\T Q_{ij} F_j;
                    lmiterm([-1, 1, n_j+1, P],-1,E_ij); 
                    matTemp2 = C_i'*S_ij + C_i'*Q_ij*F_j;
                    lmiterm([-1, 1, n_j+1, 0],matTemp2); 
                    % W_ij3 = -E_{ji}^\T P_j + S_{ji}^\T C_j + F_i^\T Q_{ji}^\T C_j;
                    matTemp3 = -E_ji'*P_j + S_ji'*C_j + F_i'*Q_ji'*C_j;
                    lmiterm([-1, 2, n_j, 0],matTemp3);
                    % W_ij4 = F_i^\T Q_{ij} F_j + (F_i^\T S_{ij} + S_{ji}^\T F_j) + R_{ij};
                    matTemp4 = F_i'*Q_ij*F_j + (F_i'*S_ij+S_ji'*F_j) + R_ij;
                    lmiterm([-1, 2, n_j+1, 0],matTemp4);
                    
                    % Filling W_ij^T = [W_ij1', W_ij3'; W_ij2', W_ij4']
                    % W_ij1' = - A_{ij}^T P_i + (-A_{ji}^T P_j + C_i^\T Q_{ij} C_j)^T;
                    lmiterm([-1, n_j, 1, P],-A_ij',1);
                    lmiterm([-1, n_j, 1, 0],matTemp1');
                    % W_ij3' = (-E_{ji}^\T P_j + S_{ji}^\T C_j + F_i^\T Q_{ji}^\T C_j)^T;
                    lmiterm([-1, n_j, 2, 0],matTemp3');
                    % W_ij2' = -E_{ij}^T P_i  + (C_i^\T S_{ij} + C_i^\T Q_{ij} F_j)^T;
                    lmiterm([-1, n_j+1, 1, P],-E_ij',1); 
                    lmiterm([-1, n_j+1, 1, 0],matTemp2'); 
                    % W_ij4' = (F_i^\T Q_{ij} F_j + (F_i^\T S_{ij} + S_{ji}^\T F_j) + R_{ij})^T
                    lmiterm([-1, n_j+1, 2, 0],matTemp4');
                    
                end

                % P term
                lmiterm([-2, 1, 1, P],1,1); % defines -P<0
                
                % Completing the LMI
                lmisys = getlmis;
                [tmin,Psol] = feasp(lmisys); % Solve the LMI system
                isFeasible = tmin <= 0; % strictly feasible if this is satisfied
                P_i = dec2mat(lmisys,Psol,P); % This needs to be stored
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
                    S_ij = obj.dataToBeDistributed.S{jInd};
                    S_ji = subsystems(jInd).dataToBeDistributed.S{iInd};
                    R_ij = obj.dataToBeDistributed.R{jInd};
                    
                    
                    % Filling W_ij = [W_ij1, W_ij2; W_ij3, W_ij4]
                    % W_ij1 = -A_{ji}^T P_j - P_iA_{ij} + C_i^\T Q_{ij} C_j;
                    W_ij1 = -P_i*A_ij -A_ji'*P_j + C_i'*Q_ij*C_j;
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
                % W_ii1 = -A_{ii}^T P_i - P_iA_{ii} + C_i^\T Q_{ii} C_i
                W_ii1 = -A_ii'*P_i - P_i*A_ii + C_i'*Q_ii*C_i;
                % W_ii2 = -P_i E_{ii} + C_i^\T S_{ii} + C_i^\T Q_{ii} F_i;
                W_ii2 = -P_i*E_ii + C_i'*S_ii + C_i'*Q_ii*F_i;
                % W_ii3 = -E_{ii}^\T P_i + S_{ii}^\T C_i + F_i^\T Q_{ii}^\T C_i;
                W_ii3 = W_ii2';
                % W_ii4 = F_i^\T Q_{ii} F_i + (F_i^\T S_{ii} + S_{ii}^\T F_i) + R_{ii};
                W_ii4 = F_i'*Q_ii*F_i + (F_i'*S_ii+S_ii'*F_i) + R_ii;
                
                W_ii = [W_ii1, W_ii2; W_ii3, W_ii4];
                
                tildeW_ii = W_ii - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                tildeW_i = [tildeW_i, tildeW_ii];

                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                if ~isFeasible
                    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            end
          
        end
        
        
        %% Local stabilizing control design LMI approach
        function [isStabilizable,K_ii,K_ijVals,K_jiVals] = designLocalStabilizingSFBControllers(obj,previousSubsystems, subsystems)
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Stabilizing at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            % W_ij = -A_{ji}^T P_j - P_iA_{ij} - K_{ji}^T B_j^T P_j - P_iB_iK_{ij}
            % K_ij is p_i X n_j
            
            A_ii = obj.A{iInd};
            B_i = obj.B{iInd};
            dim_n = obj.dim_n;
            dim_p = obj.dim_p; 
            
            threshold = 0;
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii > 0  <==>  Q' W_ii Q > 0 with Q = Q' = P_i^{-1} > 0 
                % W_ii = - A_{ii}^T P_i - P_iA_{ii} - K_{ii}^T B_i^T P_i - P_iB_iK_{ii}
                % Transformed to: M_ii = -Q_iA_{ii}^T - A_{ii}Q_i - L_i^T B_i^T - B_i L_i with L_i = K_{ii}Q_i
                % we need M_ii>0 and Q_i > 0 and  L_i - free (then P_i = Q_i^{-1}, K_ii = L_i P_i)
                
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                Q = lmivar(1,[dim_n, 1]);   % Q is a variable, 1: square symmetric, size(A,1) gives the size and 1 gives that Q is a full matrix 
                L = lmivar(2,[dim_p,dim_n])      % L is a variable, rectangular matrix
                
                % M_ii = -Q_iA_{ii}^T - A_{ii}Q_i - L_i^T B_i^T - B_i L_i
                lmiterm([-1, 1, 1, Q],-A_ii,1,'s'); % defines -(-Q_iA_{ii}^T - A_{ii}Q_i)<0; here "s" flag is used to get the symmetric part of PA.
                lmiterm([-1, 1, 1, L],-B_i,1,'s'); % defines -(- L_i^T B_i^T - B_i L_i)<0;
                lmiterm([-2, 1, 1, Q],1,1); % defines -Q<0
                lmisys = getlmis;
                
                [tmin,sol] = feasp(lmisys); % Solve the LMI system
                isStabilizable = tmin < threshold; % strictly feasible if this is satisfied
                Q_i = dec2mat(lmisys,sol,Q); 
                L_i = dec2mat(lmisys,sol,L); 
                P_i = eye(dim_n)/Q_i; % This needs to be stored
                K_ii = L_i*P_i; % This needs to be stored
                K_jiVals = [];
                K_ijVals = [];
                tildeW_i = -A_ii'*P_i - P_i*A_ii - K_ii'*B_i'*P_i - P_i*B_i*K_ii; % Note that here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                % LMI End
                
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                obj.controllerGains.localSFBStabLMI{iInd} = K_ii;
                
                if ~isStabilizable
                    disp(['Not stabilizable at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            
            else
                % This subsystem has to talk with all the previosSubsystems
                % tildeW_ii > 0 iff [W_ii, W_i; W_i', M_i] > 0 is required where 
                % M_i = inv((scriptD_i*scriptA_i)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i)^{-1}')
                % otherwise : M_i = scriptA_i'*scriptD_i'*scriptA_i 
                
                % Introduce L_ij = P_i B_i K_ij, L_ji = K_ji' as variables apart from variables P_i.
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                P = lmivar(1,[dim_n, 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    % L_ij = lmivar(2,[dim_n, dim_n])     % For L_ij = P_i B_i K_ij
                    % L_ji = lmivar(2,[dim_n, dim_p]) % For L_ji = K_ji^T
                    eval(['L_i',num2str(j),'= lmivar(2,[dim_n, dim_n]);'])
                    eval(['L_',num2str(j),'i= lmivar(2,[dim_n, dim_p]);'])
                end
                L = lmivar(2,[dim_n,dim_n]); % For L_i = P_i B_i K_ii  
                
%                 W_ii term: W_ii = - A_{ii}^T P_i - P_iA_{ii} - K_{ii}^T B_i^T P_i - P_iB_iK_{ii}
                % W_ii = - A_{ii}^T P_i - P_iA_{ii} - L' - L
                lmiterm([-1, 1, 1, P],-1,A_ii,'s'); % defines -(-PA-A'P)<0; here "s" flag is used to get the symmetric part of PA.
                lmiterm([-1, 1, 1, L],-1,1,'s'); % defines -(-PA-A'P)<0; here "s" flag is used to get the symmetric part of PA.
                
                % M_i term
                blockSize = dim_n; 
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
                disp(['Data at ',num2str(iInd),' after ',num2str(previousSubsystems),'.'])
                scriptD_i
                scriptA_i
                M1_i = inv(scriptD_i*scriptA_i);
                % M_i = inv(M1_i*scriptD_i*M1_i'); % THis fills (i-1)x(i-1) blocks in the LMI
                M_i = scriptA_i'*scriptD_i*scriptA_i;
                if issymmetric(scriptD_i) & issymmetric(scriptA_i) & ~issymmetric(M_i)
                    tf = norm(M_i-M_i.',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M_i = 0.5*(M_i + M_i');
                end
                
                % Filling M_i
                for x=1:1:length(previousSubsystems)
                    m_x = (x-1)*blockSize;  % exact index (for data extraction)
                    n_x = 1+(x-1)+1;        % block index (for LMI blocks)
                    for y=1:1:length(previousSubsystems)
                        m_y = (y-1)*blockSize;      % exact index (for data extraction)
                        n_y = 1+(y-1)+1;            % block index (for, LMI blocks)
                        M_i_xy = M_i(m_x+1:m_x+blockSize, m_y+1:m_y+blockSize);
                        lmiterm([-1, n_x, n_y, 0],M_i_xy);
                    end
                end
                
                
                % W_i term
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % W_ij term: -A_{ji}^T P_j - P_iA_{ij} - K_{ji}^T B_j^T P_j - P_iB_iK_{ij}
                    % W_ij = -A_{ji}^T P_j - P_iA_{ij} - L_ji B_j^T P_j - L_ij
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};
                    B_j  = subsystems(jInd).B{jInd};
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    n_j = 1+(j-1)+1; % block index (for LMI blocks)
                    
                    % Filling W_ij = -A_{ji}^T P_j - P_iA_{ij} - L_ji B_j^T P_j - L_ij
                    lmiterm([-1, 1, n_j, P],-1,A_ij);
                    matTemp1 = -A_ji'*P_j;
                    lmiterm([-1, 1, n_j, 0],matTemp1);
                    matTemp2 = B_j'*P_j;
                    % lmiterm([-1, 1, n_j, L_ji],-1,matTemp2);
                    % lmiterm([-1, 1, n_j, L_ij],-1,1);
                    eval(['lmiterm([-1, 1, n_j, L_',num2str(j),'i],-1,matTemp2);']);
                    eval(['lmiterm([-1, 1, n_j, L_i',num2str(j),'],-1,1);']);
                    
                    % Filling W_ij^T = -P_j A_{ji} - A_{ij}^T P_i - P_j B_j L_ji^T - L_ij^T
                    lmiterm([-1, n_j, 1, P],-A_ij',1);
                    lmiterm([-1, n_j, 1, 0],matTemp1');
                    % lmiterm([-1, n_j, 1, -L_ji],-matTemp2',1);
                    % lmiterm([-1, n_j, 1, -L_ij],-1,1);
                    eval(['lmiterm([-1, n_j, 1, -L_',num2str(j),'i],-transpose(matTemp2),1);']);
                    eval(['lmiterm([-1, n_j, 1, -L_i',num2str(j),'],-1,1);']);
                end
                
                % P term
                lmiterm([-2, 1, 1, P],1,1); % defines -P<0
                
                % Completing the LMI
                lmisys = getlmis;
                [tmin,sol] = feasp(lmisys); % Solve the LMI system
                isStabilizable = tmin < threshold; % strictly feasible if this is satisfied
                P_i = dec2mat(lmisys,sol,P); % This needs to be stored
%                 B_i
                L_i = dec2mat(lmisys,sol,L); % L_i = P_i B_i K_ii
                K_ii = pinv(B_i)*inv(P_i)*L_i;  % same as (B_i'*B_i)\B_i'*inv(P_i)*L_i
                obj.controllerGains.localSFBStabLMI{iInd} = K_ii;
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    % L_ij = P_i B_i K_ij
                    % L_ij = dec2mat(lmisys,sol,L_ij); 
                    eval(['L_ij = dec2mat(lmisys,sol,L_i',num2str(j),');']);
                    K_ij = pinv(B_i)*inv(P_i)*L_ij;
                    K_ijVals{jInd} = K_ij;
                    L_ijVals{jInd} = L_ij;
                    obj.controllerGains.localSFBStabLMI{jInd} = K_ij;
                    
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
                    B_j  = subsystems(jInd).B{jInd};
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    K_ij = K_ijVals{jInd}; % same as obj.localStabilizingSFBControllerGains{jInd};
                    K_ji = K_jiVals{jInd};
                    L_ij = L_ijVals{jInd};
                    L_ji = L_jiVals{jInd};
                    % Filling W_ij = -A_{ji}^T P_j - P_iA_{ij} - L_ji B_j^T P_j - L_ij
                    % W_ij = -A_{ji}^T P_j - P_iA_{ij} - K_{ji}^T B_j^T P_j - P_iB_iK_{ij}
%                     W_ij = -P_i*A_ij - A_ji'*P_j - K_ji'*B_j'*P_j - P_i*B_i*K_ij
                    W_ij = -P_i*A_ij - A_ji'*P_j - L_ji*B_j'*P_j - L_ij;
%                     W_ijE-W_ij
%                     norm(W_ijE-W_ij)
                    W_i = [W_i, W_ij];
                end
                tildeW_i = W_i*M1_i;
                
                % W_ii term: - A_{ii}^T P_i - P_iA_{ii} - K_{ii}^T B_i^T P_i - P_iB_iK_{ii}
%                 W_ii = -A_ii'*P_i - P_i*A_ii - K_ii'*B_i'*P_i - P_i*B_i*K_ii
                W_ii = -A_ii'*P_i - P_i*A_ii - L_i' - L_i;
%                 W_iiE-W_ii
%                 norm(W_iiE-W_ii)
                
                tildeW_ii = W_ii - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                tildeW_i = [tildeW_i, tildeW_ii];

                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                if ~isStabilizable
                    disp(['Not stabilizable at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            end
        end
        
        
        
        %% The alternative approach to local stabilizing controller design using LMIs
        function [isStabilizable,K_ii,K_ijVals,K_jiVals] = designLocalStabilizingSFBControllers2(obj,previousSubsystems, subsystems)
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Stabilizing at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            % K_ij is p_i X n_j
            
            A_ii = obj.A{iInd};
            B_i = obj.B{iInd};
            dim_n = obj.dim_n;
            dim_p = obj.dim_p; 
            
            threshold = 0;
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii > 0
                % W_ii = - A_{ii}^T P_i - P_iA_{ii} - K_{ii}^T B_i^T P_i - P_iB_iK_{ii}
                % Transformed to: M_ii = Q_i M_1^T + M_1 Q_i + L_i^T M_2^T + M_2 L_i with L_i = K_{ii}Q_i
                % we need M_ii>0 and Q_i > 0 and  L_i - free (then P_i = Q_i^{-1}, K_ii = L_i P_i)
                
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                Q = lmivar(1,[dim_n, 1]);   % Q is a variable, 1: square symmetric, size(A,1) gives the size and 1 gives that Q is a full matrix 
                L = lmivar(2,[dim_p,dim_n]);      % L is a variable, rectangular matrix
                
                % M_ii = Q_i M_1^T + M_1 Q_i + L_i^T M_2^T + M_2 L_i
                lmiterm([-1, 1, 1, Q],-A_ii,1,'s'); % defines -(-Q_iA_{ii}^T - A_{ii}Q_i)<0; here "s" flag is used to get the symmetric part of PA.
                lmiterm([-1, 1, 1, L],-B_i,1,'s'); % defines -(- L_i^T B_i^T - B_i L_i)<0;
                lmiterm([-2, 1, 1, Q],1,1); % defines -Q<0
                lmisys = getlmis;
                
                [tmin,sol] = feasp(lmisys); % Solve the LMI system
                isStabilizable = tmin < threshold; % strictly feasible if this is satisfied
                Q_i = dec2mat(lmisys,sol,Q); 
                L_i = dec2mat(lmisys,sol,L); 
                P_i = eye(dim_n)/Q_i; % This needs to be stored
                K_ii = L_i*P_i; % This needs to be stored
                K_jiVals = [];
                K_ijVals = [];
                tildeW_i = -A_ii'*P_i - P_i*A_ii - K_ii'*B_i'*P_i - P_i*B_i*K_ii; % Note that here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                % LMI End
                
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                obj.controllerGains.localSFBStabLMI{iInd} = K_ii;
                
                if ~isStabilizable
                    disp(['Not stabilizable at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            
            else
                % This subsystem has to talk with all the previosSubsystems
                % tildeW_ii > 0 iff [W_ii, W_i; W_i', M_i] > 0 is required where                 
                
                % LMI Start
                setlmis([]);  % To initialize the LMI description
                
                % % Introduce LMI variables: L_i and Q_i and {L_ji:j=1,2,..,i-1} .
                
                % Defining Q_i = diag(inv(P_i),I,I,...)
                [I_nn,~,sI_nn] = lmivar(1,[dim_n, 0]); % diagonal
                [Qi1,~,sQi1] = lmivar(1,[dim_n, 1]);   % symmetric
                sQi = [sQi1];
                for k = 1:1:(i-1)
                    sQi = [sQi, zeros(k*dim_n,dim_n); zeros(dim_n,k*dim_n), sI_nn]; %
                end
                [Qi,~,sQi] = lmivar(3,sQi);
                sQi
                
                % Defining L_i = diag(L_ii,L_i1,L_i2,...,L_i,i-1) with L_ii = K_iiP_i^{-1}, L_i1 = K_i1, L_i2 = K_i2,...
                [Lij,~,sLij] = lmivar(2,[dim_p,dim_n]); % rectangular
                sLi = [sLij];
                for j = 1:1:(i-1)
                    % [Lij,~,sLij] = lmivar(2,[dim_p dim_n]);
                    eval(['[Li',num2str(j),',~,sLij] = lmivar(2,[dim_p dim_n]);']);
                    sLi = [sLi, zeros(j*dim_p,dim_n); zeros(dim_p,j*dim_n), sLij];
                end
                [Li,~,sLi] = lmivar(3,sLi);
                sLi
                
                % Defining L_ji = diag(L_ji1,I,I,...) with L_ji1 = K_ji*inv(P_i)
                for j = 1:1:(i-1)
                    [Lji1,~,sLji1] = lmivar(2,[dim_p dim_n]); % rectangular
                    sLji = [sLji1, zeros(dim_p,(i-1)*dim_n)];
                    for k = 1:1:(i-1)
                        sLji = [sLji;[zeros(dim_n,k*dim_n), sI_nn, zeros(dim_n,(i-1-k)*dim_n)]];
                    end
                    % [Lji,~,sLji] = lmivar(3,sLji);
                    eval(['[L',num2str(j),'i,~,sL',num2str(j),'i] = lmivar(3,sLji);'])
                    eval(['sL',num2str(j),'i'])
                end
                
                
                % % Introduce constants: M1,M2, and M_i, M3 and {M4j:j=1,2,...,i-1}
                
                % Defining M1 = [-A_ii,-A_i1,-A_i2,...,-A_i,i-1; 0,0,0,...,0; 0,0,0,...,0;...]
                Z1 = zeros(dim_n*(i-1), dim_n);
                M1 = [-A_ii; Z1];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    
                    A_ij = obj.A{jInd};
                    M1 = [M1, [-A_ij; Z1]];
                end
                M1
                
                % Defining M2 = [-B_i,-B_i,...,-B_i; 0,0,0,...,0; 0,0,0,...,0;...]
                Z2 = zeros(dim_n*(i-1), dim_p);
                M2 = [-B_i; Z2];
                for j = 1:1:length(previousSubsystems)
                    M2 = [M2, [-B_i; Z2]];
                end
                M2
                
                
                % Defining M based on the information recieved (following the same old approach)
                blockSize = dim_n; 
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
                disp(['Data at ',num2str(iInd),' after ',num2str(previousSubsystems),'.'])
                % scriptD_i
                % scriptA_i
                M1_i = inv(scriptD_i*scriptA_i);
                % M_i = inv(M1_i*scriptD_i*M1_i'); % THis fills (i-1)x(i-1) blocks in the LMI
                M_i = scriptA_i'*scriptD_i*scriptA_i;
                if issymmetric(scriptD_i) & issymmetric(scriptA_i) & ~issymmetric(M_i)
                    tf = norm(M_i-M_i',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M_i = 0.5*(M_i + M_i');
                end
                
                
                % Defining M3 = inv([Z3,M31;M31',M/i])
                Z3 = zeros(dim_n,dim_n);
                M31 = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    A_ji = subsystems(jInd).A{iInd};
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    
                    M31 = [M31, -A_ji'*P_j];
                end
                M3inv = [Z3,M31;M31',M_i/i]
                if all(M31(:)==0) % if M31 is zero, remove M_i/i from here and take M3 = 0
                    M3 = 0*M3inv % all zeros
                    M_iBlock = M_i/(i-1);
                else
                    M3 = inv(M3inv)
                    M_iBlock = M_i/i;
                end
                if ~issymmetric(M3)
                    tf = norm(M3-M3',inf);
                    disp(['Symmetry Error M3 !!! Magnitude:',num2str(tf)]);
                    M3 = 0.5*(M3 + M3');
                end               
                
                % Defining M4j for j=1,2,...,i-1 with M4j = inv([Z4,M4j1;M4j1',M/i])
                Z4 = zeros(dim_p,dim_p);
                for j = 1:1:length(previousSubsystems)
                    % Defining M4j = inv([Z4,M4j1;M4j1',M_i/i])
                    jInd = previousSubsystems(j);
                    B_j = subsystems(jInd).B{jInd}
                    P_j = subsystems(jInd).dataToBeDistributed.P
                    
                    Z4j = zeros(dim_p,(j-1)*dim_n);
                    z4j = zeros(dim_p,(i-1-j)*dim_n);
                    M4j1 = [Z4j,-B_j'*P_j,z4j];
                
                    M4jinv = [Z4, M4j1;M4j1',M_iBlock]
                    %det(M4j1*inv(M_iBlock)*M4j1')
                    if abs(det(M4jinv))<0.000000001
                        M4j = 0*M4jinv
                    else
                        M4j = inv(M4jinv);
                        if ~issymmetric(M4j)
                            tf = norm(M4j-M4j',inf);
                            disp(['Symmetry Error M4j !!! Magnitude:',num2str(tf)]);
                            M4j = 0.5*(M4j + M4j');
                        end
                    end
                    eval(['M4',num2str(j),' = M4j;'])
                    eval(['M4',num2str(j)])
                    % M4j = inv([Z4, M4j1;M4j1',M_iBlock]);
                    % eval(['M4',num2str(j),' = inv([Z4, M4j1;transpose(M4j1),M_iBlock])'])
                end
                
                
                % % Adding LMI Terms 
                lmiterm([-1, 1, 1, Qi],M1,1,'s');
                lmiterm([-1, 1, 1, Li],M2,1,'s');
                lmiterm([-1, 1, 2, -Qi],1,1);
                lmiterm([-1, 2, 1, Qi],1,1);
                lmiterm([-1, 2, 2, 0],M3);
                Z5 = zeros(size(M3,1),size(M41,2));
                for j = 1:1:length(previousSubsystems)
                    % lmiterm([-1, 1, 2+j, -Lji],1,1);
                    eval(['lmiterm([-1, 1, 2+j, -L',num2str(j),'i],1,1);']);
                    lmiterm([-1, 2, 2+j, 0],Z5);
                    % lmiterm([-1, 2+j, 1, Lji],1,1);
                    eval(['lmiterm([-1, 2+j, 1, L',num2str(j),'i],1,1);']);
                    lmiterm([-1, 2+j, 2, 0],Z5');
                end                
                
                if length(previousSubsystems)>=2
                    Z6 = zeros(size(M41,1),size(M42,2));
                end
                for j = 1:1:length(previousSubsystems)
                    for k = 1:1:length(previousSubsystems)
                        if j==k
                            %lmiterm([-1, 2+j, 2+k, 0],M4j);
                            eval(['lmiterm([-1, 2+j, 2+k, 0],M4',num2str(j),');']);
                        elseif j<k
                            lmiterm([-1, 2+j, 2+k, 0],Z6);
                        else
                            lmiterm([-1, 2+j, 2+k, 0],Z6');
                        end
                    end
                end
                
                lmiterm([-2, 1, 1, Qi],1,1);
                
                % % Completing the LMI
                lmisys = getlmis;
                lmisys = setmvar(lmisys,I_nn,1);
                [tmin,sol] = feasp(lmisys);
                
% %                 evals = evallmi(lmisys,sol);
% %                 [lhs1,rhs1] = showlmi(evals,1) 
% %                 [lhs2,rhs2] = showlmi(evals,2)
% %                 eig1 = eig(lhs1-rhs1)
% %                 eig2 = eig(lhs2-rhs2)
% %                 eig1(eig1>0)
% %                 all(eig1<0)
% %                 all(eig2<0)
                
                isStabilizable = tmin < threshold % strictly feasible if this is satisfied
                
                % % decoding solutions
                Qisol = dec2mat(lmisys,sol,Qi) % for some reason preset variables not appear in Q (they appear as zeros)
                P_i = inv(Qisol(1:dim_n,1:dim_n))
                Lisol = dec2mat(lmisys,sol,Li)
                K_ii = Lisol(1:dim_p,1:dim_n)*P_i
                obj.controllerGains.localSFBStabLMI{iInd} = K_ii;
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    K_ij = Lisol(j*dim_p+1:(j+1)*dim_p, j*dim_n+1:(j+1)*dim_n)
                    K_ijVals{jInd} = K_ij;
                    
                    % Ljisol = dec2mat(lmisys,sol,Lji)
                    eval(['Ljisol = dec2mat(lmisys,sol,L',num2str(j),'i)']);
                    K_ji = Ljisol(1:dim_p,1:dim_n)*P_i
                    K_jiVals{jInd} = K_ji;
                end
                
                
                % % Finally, we need to compute \tilede{W}_i and \tilde{W}_{ii} for storage
                
                % W_i term
                W_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % W_ij term
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};
                    B_j  = subsystems(jInd).B{jInd};
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    K_ij = K_ijVals{jInd}; % same as obj.localStabilizingSFBControllerGains{jInd};
                    K_ji = K_jiVals{jInd};
                    % Filling W_ij = -A_{ji}^T P_j - P_iA_{ij} - L_ji B_j^T P_j - L_ij
                    % W_ij = -A_{ji}^T P_j - P_iA_{ij} - K_{ji}^T B_j^T P_j - P_iB_iK_{ij}
                    W_ij = - P_i*A_ij - A_ji'*P_j - K_ji'*B_j'*P_j - P_i*B_i*K_ij
                    
%                     W_ijE-W_ij
%                     norm(W_ijE-W_ij)
                    W_i = [W_i, W_ij];
                end
                tildeW_i = W_i*M1_i;
                
                % W_ii term: - A_{ii}^T P_i - P_iA_{ii} - K_{ii}^T B_i^T P_i - P_iB_iK_{ii}
                W_ii = -A_ii'*P_i - P_i*A_ii - K_ii'*B_i'*P_i - P_i*B_i*K_ii
%                 W_iiE-W_ii
%                 norm(W_iiE-W_ii)
                
                tildeW_ii = W_ii - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                tildeW_i = [tildeW_i, tildeW_ii];

                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                if ~isStabilizable
                    disp(['Not stabilizable at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            end
        end
        
        
        
        %% ALternative approach to Local stabilizing control design using YALMIP 
        function [isStabilizable,K_ii,K_ijVals,K_jiVals] = designLocalStabilizingSFBControllers3(obj,previousSubsystems, subsystems)
            i = length(previousSubsystems)+1;
            iInd = obj.index;
            disp(['Stabilizing at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
            % W_ij = -A_{ji}^T P_j - P_iA_{ij} - K_{ji}^T B_j^T P_j - P_iB_iK_{ij}
            % K_ij is p_i X n_j
            
            A_ii = obj.A{iInd};
            B_i = obj.B{iInd};
            dim_n = obj.dim_n;
            dim_p = obj.dim_p; 
            
            if isempty(previousSubsystems)
                % The subsystem only need to test W_ii > 0  <==>  Q' W_ii Q > 0 with Q = Q' = P_i^{-1} > 0 
                % W_ii = - A_{ii}^T P_i - P_iA_{ii} - K_{ii}^T B_i^T P_i - P_iB_iK_{ii}
                % Transformed to: M_ii = -Q_iA_{ii}^T - A_{ii}Q_i - L_i^T B_i^T - B_i L_i with L_i = K_{ii}Q_i
                % we need M_ii>0 and Q_i > 0 and  L_i - free (then P_i = Q_i^{-1}, K_ii = L_i P_i)
                
                % LMI Start
                % setlmis([]);  % To initialize the LMI description
                % Q = lmivar(1,[dim_n, 1]);   % Q is a variable, 1: square symmetric, size(A,1) gives the size and 1 gives that Q is a full matrix 
                Q = sdpvar(dim_n,dim_n);
                % L = lmivar(2,[dim_p,dim_n])      % L is a variable, rectangular matrix
                L = sdpvar(dim_p,dim_n,'full');
                
                % M_ii = -Q_iA_{ii}^T - A_{ii}Q_i - L_i^T B_i^T - B_i L_i
                M_ii = - Q*A_ii' - A_ii*Q - L'*B_i' - B_i*L;
                constraints = [Q>=0, M_ii>=0];
                sol = optimize(constraints);
                isStabilizable = sol.problem==0;
                
                Q_i = value(Q);
                L_i = value(L);
                P_i = eye(dim_n)/Q_i; % This needs to be stored
                K_ii = L_i*P_i; % This needs to be stored
                K_jiVals = [];
                K_ijVals = [];
                tildeW_i = -A_ii'*P_i - P_i*A_ii - K_ii'*B_i'*P_i - P_i*B_i*K_ii; % Note that here, \tilde{W}_ii = W_ii = \tilde{W}_i. This also needs to be stored
                % LMI End
                
                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                obj.controllerGains.localSFBStabLMI{iInd} = K_ii;
                
                if ~isStabilizable
                    disp(['Not stabilizable at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                end
            
            else
                % This subsystem has to talk with all the previosSubsystems
                % tildeW_ii > 0 iff [W_ii, W_i; W_i', M_i] > 0 is required where 
                % M_i = inv((scriptD_i*scriptA_i)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i)^{-1}')
                % otherwise : M_i = scriptA_i'*scriptD_i'*scriptA_i 
                
                % Introduce L_ij = P_i B_i K_ij, L_ji = K_ji' as variables apart from variables P_i.
                % LMI Start
                % setlmis([]);  % To initialize the LMI description
                % P = lmivar(1,[dim_n, 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 
                P_i = sdpvar(dim_n,dim_n);
                % L = lmivar(2,[dim_n,dim_n]); % For L_i = P_i B_i K_ii  
                K_ii = sdpvar(dim_p,dim_n,'full');
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    % L_ij = lmivar(2,[dim_n, dim_n])     % For L_ij = P_i B_i K_ij
                    % L_ji = lmivar(2,[dim_n, dim_p]) % For L_ji = K_ji^T
                    K_ij = sdpvar(dim_p,dim_n,'full');
                    K_ji = sdpvar(dim_p,dim_n,'full');
                    eval(['K_i',num2str(j),' = sdpvar(dim_p,dim_n);']);
                    eval(['K_',num2str(j),'i = sdpvar(dim_p,dim_n);']);
                end
                
                
                % W_ii term: W_ii = - A_{ii}^T P_i - P_iA_{ii} - K_{ii}^T B_i^T P_i - P_iB_iK_{ii}
                W_ii = -A_ii'*P_i - P_i*A_ii - K_ii'*B_i'*P_i - P_i*B_i*K_ii;
                
                
                % M_i term:
                blockSize = dim_n; 
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
                disp(['Data at ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
                scriptD_i
                scriptA_i
                M1_i = inv(scriptD_i*scriptA_i);
                % M_i = inv(M1_i*scriptD_i*M1_i'); % THis fills (i-1)x(i-1) blocks in the LMI
                M_i = scriptA_i'*scriptD_i*scriptA_i;
                if issymmetric(scriptD_i) & issymmetric(scriptA_i) & ~issymmetric(M_i)
                    tf = norm(M_i-M_i.',inf);
                    disp(['Symmetry Error !!! Magnitude:',num2str(tf)]);
                    % M_i
                    M_i = 0.5*(M_i + M_i');
                end
                   
                
                % W_i term:
                W_i = [];
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);

                    % W_ij term: -A_{ji}^T P_j - P_iA_{ij} - K_{ji}^T B_j^T P_j - P_iB_iK_{ij}
                    A_ij = obj.A{jInd};
                    A_ji = subsystems(jInd).A{iInd};
                    B_j  = subsystems(jInd).B{jInd};
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    n_j = 1+(j-1)+1; % block index (for LMI blocks)
                    
                    % Filling W_ij = -A_{ji}^T P_j - P_iA_{ij} - K_{ji}^T B_j^T P_j - P_iB_iK_{ij}
                    W_ij = - A_ji'*P_j - P_i*A_ij;                    
                    W_ij = W_ij + eval(['-transpose(K_',num2str(j),'i)*transpose(B_j)*P_j']);
                    W_ij = W_ij + eval(['P_i*B_i*K_i',num2str(j)]);
                    
                    W_i = [W_i, W_ij];                    
                end
                
                
                % Constaints
                constraints = [P_i>=0, [W_ii,W_i;W_i',M_i]>=0];
                
                
                % Completing the LMI
                sol = optimize(constraints);
                isStabilizable = sol.problem==0;
                
                P_i = value(P_i);
                K_ii = value(K_ii);
                obj.controllerGains.localSFBStabLMI{iInd} = K_ii;
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    K_ij = value(eval(['K_i',num2str(j)]));
                    K_ji = value(eval(['K_',num2str(j),'i']));
                    K_ijVals{jInd} = K_ij;
                    obj.controllerGains.localSFBStabLMI{jInd} = K_ij;
                    K_jiVals{jInd} = K_ji; % these values will be loaded outside the function
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
                    B_j  = subsystems(jInd).B{jInd};
                    P_j = subsystems(jInd).dataToBeDistributed.P;
                    K_ij = K_ijVals{jInd}; % same as obj.localStabilizingSFBControllerGains{jInd};
                    K_ji = K_jiVals{jInd};
                    % Filling W_ij = -A_{ji}^T P_j - P_iA_{ij} - K_{ji}^T B_j^T P_j - P_iB_iK_{ij}
                    W_ij = -P_i*A_ij - A_ji'*P_j - K_ji'*B_j'*P_j - P_i*B_i*K_ij;

                    W_i = [W_i, W_ij];
                end
                tildeW_i = W_i*M1_i;
                
                % W_ii term: - A_{ii}^T P_i - P_iA_{ii} - K_{ii}^T B_i^T P_i - P_iB_iK_{ii}
                W_ii = -A_ii'*P_i - P_i*A_ii - K_ii'*B_i'*P_i - P_i*B_i*K_ii;

                
                tildeW_ii = W_ii - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
                tildeW_i = [tildeW_i, tildeW_ii];

                obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
                obj.dataToBeDistributed.P = P_i; % Storing
                if ~isStabilizable
                    disp(['Not stabilizable at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
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

