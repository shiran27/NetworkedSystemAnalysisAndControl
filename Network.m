classdef Network < handle
    % This is the class for networks (networked systems)
        
    properties
        index
        subsystems = []         % list of subsystem objects
        edges = []              % list of edge objects
        distanceMatrix
        
        neighbors = {}          % \bar{\mathcal{E}}_i
        outNeighbors = {}       % \bar{\mathcal{F}}_i
        exclusiveNeighbors = []
        exclusiveOutNeighbors = []
        
        % With respect to a given indexing scheme
        leastNeighbors = []     % \Gamma_i 's
        leastOutNeighbors = []  % \Delta_i 's
        
        % networkMatrices
        networkMatrices
        
        % test network matrix
        testNetworkMatrix
    end
    
    methods
        
        % Constructor
        function obj = Network(index)
            % Construct an instance of this class
            obj.index = index;
        end
        
        
        % function to load a random networked system
        function obj = loadARandomNetwork(obj,numOfSubsystems,dimentionOfSpace,sizeOfSpace,communicationRadius)
            
            % Generate subsystems' list
            subsystemLocations = sizeOfSpace*rand(numOfSubsystems,dimentionOfSpace);
            
            % Initializing subsystems: 
            % Subsystem-spaecific parameters will be loaded inside
            obj.subsystems = [];
            for i = 1:1:numOfSubsystems
                newSubsystem = Subsystem(i, subsystemLocations(i,:));
                obj.subsystems = [obj.subsystems, newSubsystem];
            end
            
        
            % Generate edges' list
            edgeCount = 1;
            obj.edges = [];
            for i = 1:1:numOfSubsystems
                for j = 1:1:numOfSubsystems
                    if i~=j
                        newEdge = Edge(edgeCount,[obj.subsystems(i).index, obj.subsystems(j).index],[obj.subsystems(i).location; obj.subsystems(j).location]);
                        if rand(1)<0.5 & norm(obj.subsystems(i).location - obj.subsystems(j).location, 2) < communicationRadius
                            newEdge.enabled = true;
                        else
                            newEdge.enabled = false;
                        end
                        obj.edges = [obj.edges, newEdge];
                        edgeCount = edgeCount + 1;
                    end
                end
            end
            
            % loading distances matrix and neighbors
            obj = obj.loadDistanceMatrix();
            obj.loadNeighbors();
            
            
            % Loading subsystem parameters
            for i = 1:1:numOfSubsystems
%                 obj.subsystems(i).loadParameters(obj.subsystems);
                obj.subsystems(i).loadStableParameters(obj.subsystems);
                obj.subsystems(i).designLocalSFBLQRControllerGains(); % find local state feedback LQR controller gains
                
            end
            
            
        end
        
        % function to load a custom networked system
        function obj = loadTheCustomNetwork(obj)
            
            % Subsystem locations (for plotting purposes only) Nx2 matrix
            subsystemLocations = [0.5,0.5; 0.2, 0.2; 0.8, 0.2; 0.8,0.8; 0.2,0.8];
            numOfSubsystems = size(subsystemLocations,1);
            
            % Initializing subsystems: 
            % Subsystem-spaecific parameters will be loaded inside
            obj.subsystems = [];
            for i = 1:1:numOfSubsystems
                newSubsystem = Subsystem(i, subsystemLocations(i,:));
                obj.subsystems = [obj.subsystems, newSubsystem];
            end
            
            % Inter-subsystem connections: defined by a list of edges Nx2
            edgeList = [1,2;2,1;1,3;3,1;1,4;4,1;1,5;5,1;2,3;3,4;4,5;5,2];
        
            % Generate edges' list
            obj.edges = [];
            edgeCount = 1;
            for i = 1:1:numOfSubsystems
                for j = 1:1:numOfSubsystems
                    if i~=j
                        newEdge = Edge(edgeCount,[obj.subsystems(i).index, obj.subsystems(j).index],[obj.subsystems(i).location; obj.subsystems(j).location]);
                        
                        foundEdge = false;
                        for k = 1:1:size(edgeList,1)
                            if i==edgeList(k,1) &&  j==edgeList(k,2)
                                foundEdge = true;
                                break
                            end
                        end
                        
                        if foundEdge
                            newEdge.enabled = true;
                        else
                            newEdge.enabled = false;
                        end
                        
                        obj.edges = [obj.edges, newEdge];
                        edgeCount = edgeCount + 1;
                    end
                end
            end
            
            % loading distances matrix and neighbors
            obj = obj.loadDistanceMatrix();
            obj.loadNeighbors();
            
            
            % Loading subsystem parameters
            for i = 1:1:numOfSubsystems
%                 obj.subsystems(i).loadParameters(obj.subsystems);
                obj.subsystems(i).loadStableParameters(obj.subsystems);
                obj.subsystems(i).designLocalSFBLQRControllerGains(); % find local state feedback LQR controller gains     
            end
            
            
        end
        
        
        % To get the distance matrix
        function obj = loadDistanceMatrix(obj) 
            numOfSubsystems = length(obj.subsystems);
            
            dmat = ones(numOfSubsystems)-eye(numOfSubsystems);
            dmat(dmat==1)=1000;
            for i = 1:1:length(obj.edges)
                iInd = obj.getSubsystem(obj.edges(i).subsystems(1),'id');
                jInd = obj.getSubsystem(obj.edges(i).subsystems(2),'id');
                distVal = obj.edges(i).getLength();
                dmat(iInd,jInd) = distVal;             
            end
            obj.distanceMatrix = dmat;
        end
        
        % get a subsystem or its id
        function result = getSubsystem(obj, subsystemIndex, requirement)
            count1 = 1;
            for subsystem = obj.subsystems % search throguh all the subsystems
                if subsystem.index == subsystemIndex 
                   if isequal(requirement,'id')
                       result = count1; break;
                   elseif isequal(requirement,'subsystem')
                       result = subsystem; break;
                   end
                end
                count1 = count1 + 1;
            end
        end
        
        % load neighbor sets based on the distance matrix
        function output = loadNeighbors(obj)
            
            for currentSubsystemId = 1:1:length(obj.subsystems) 
                neighborSet = []; % conventional in-neighbors
                outNeighborSet = []; % out-neigbors
                for subsystemId = 1:1:length(obj.subsystems)
                    % distance will be large if not connected
                    if obj.distanceMatrix(subsystemId,currentSubsystemId) < 1.5  
                        neighborSet = [neighborSet, subsystemId];
                    end
                    
                    if obj.distanceMatrix(currentSubsystemId,subsystemId) < 1.5
                        outNeighborSet = [outNeighborSet, subsystemId];
                    end
                    
                end
                obj.neighbors{currentSubsystemId} = neighborSet;
                obj.subsystems(currentSubsystemId).neighbors = neighborSet;
                obj.exclusiveNeighbors{currentSubsystemId} = neighborSet(neighborSet~=currentSubsystemId);
                obj.subsystems(currentSubsystemId).exclusiveNeighbors = neighborSet(neighborSet~=currentSubsystemId);
                
                obj.outNeighbors{currentSubsystemId} = outNeighborSet;
                obj.subsystems(currentSubsystemId).outNeighbors = outNeighborSet;
                obj.exclusiveOutNeighbors{currentSubsystemId} = outNeighborSet(outNeighborSet~=currentSubsystemId);
                obj.subsystems(currentSubsystemId).exclusiveOutNeighbors = outNeighborSet(outNeighborSet~=currentSubsystemId);
                
                % least neighbors depend on the indexing we used! indexing = original indexes is assumed in the following
%                 leastNeighbor = min(neighborSet(neighborSet~=currentSubsystemId));
%                 obj.leastNeighbors(currentSubsystemId) = leastNeighbor;
%                 obj.subsystems(currentSubsystemId).leastNeighbor = leastNeighbor;
%                 
%                 leastOutNeighbor = min(outNeighborSet(outNeighborSet~=currentSubsystemId));
%                 obj.leastOutNeighbors(currentSubsystemId) = leastOutNeighbor;
%                 obj.subsystems(currentSubsystemId).leastOutNeighbor = leastOutNeighbor;
            end
        
        end
        
        
        % create a random symmetric network matrix
        function networkMatrix = getARandomNetworkMatrix(obj,blockSize)
            networkMatrix = [];
            for i = 1:1:length(obj.subsystems)
                netRow = [];
                for j = 1:1:length(obj.subsystems)
                    if any(obj.neighbors{i}==j)
                        if i==j & rand(1)>0.3
                            netRow = [netRow, 10*eye(blockSize)]; % to make pdf
                        else
                            netRow = [netRow, 5*(rand(blockSize)-0.5)];
                        end
                    else
                        netRow = [netRow, zeros(blockSize)];
                    end
                end
                networkMatrix = [networkMatrix; netRow];
            end
            networkMatrix = networkMatrix + networkMatrix'; % to make this a symmetric matrix!
        end
        
        
        % Check the positivedefiniteness of a network matrix
        function isPositiveDefinite = checkPositiveDefiniteness(obj,testNetworkMatrix,indexing)
            
            blockSize = size(testNetworkMatrix,1)/length(obj.subsystems);
            
            % default indexing
            if isempty(indexing)
                indexing = [1:1:length(obj.subsystems)];
            end
            
            % distribute testNetworkMatrix information among the agents
            for i = 1:1:length(obj.subsystems)
                for j = 1:1:length(obj.subsystems)
                    % Assigning the i,j the block
                    obj.subsystems(i).testMatrix{j} = testNetworkMatrix(blockSize*(i-1)+1:blockSize*i,blockSize*(j-1)+1:blockSize*j);
                end
            end
            
            % Execute the distributed startegy
            for i = 1:1:length(indexing)
                iInd = indexing(i);
                previousSubsystems = indexing(1:i-1);                
                isPositiveDefinite = obj.subsystems(iInd).checkPositiveDefiniteness(previousSubsystems, obj.subsystems);
                if ~isPositiveDefinite
                    output = false;
                    break
                else
                    output = true;
                end
            end
            
        end
        
        
        
        % Draw the network in matlab figure
        function outputArg = drawNetwork(obj, figNum)
            figure(figNum)
            
            grid on
            axis equal
            sizeOfSpace = 1;
            axis([0,sizeOfSpace,0,sizeOfSpace]);
            
            % draw edges
            for i = 1:1:length(obj.edges)
                obj.edges(i).drawEdge();
            end
            
            % draw subsystems
            for i = 1:1:length(obj.subsystems)
                obj.subsystems(i).drawSubsystem();
            end
            
            xlabel('X-Location')
            ylabel('Y-Location')
            axis equal
            axis([0,sizeOfSpace,0,sizeOfSpace]);
            
            % Printing Legend
            posX = 0.9;
            posY = 0.9;
            thickness = 0.05;
            height = 0.01;

            rx = rectangle('Position',[posX, posY+6*height, thickness, height],...
                'FaceColor',[0.8 0 0 0.5],'EdgeColor','k','LineWidth',0.01);
            tx = text(posX-3*height,posY+6*height+height/2,'x_i','Color','k','FontSize',10);

            ry = rectangle('Position',[posX, posY+4*height, thickness, height],...
                'FaceColor',[0 0.8 0 0.5],'EdgeColor','k','LineWidth',0.01);
            ty = text(posX-3*height,posY+4*height+height/2,'y_i','Color','k','FontSize',10);

            ru = rectangle('Position',[posX, posY+2*height, thickness, height],...
                'FaceColor',[0 0 0.8 0.5],'EdgeColor','k','LineWidth',0.01);
            ty = text(posX-3*height,posY+2*height+height/2,'u_i','Color','k','FontSize',10);

            rw = rectangle('Position',[posX, posY, thickness, height],...
                'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k','LineWidth',0.01);
            ty = text(posX-3*height,posY+height/2,'w_i','Color','k','FontSize',10);
        end
        
        
        % Generate the complete system.
        function [A,B,C,D,E,F,x] = getNetworkMatrices(obj)
            A = []; B = []; C = []; D = []; E = []; F = []; x = [];
            for i=1:1:length(obj.subsystems)
                A_i = []; B_i = []; C_i = []; D_i = []; E_i = []; F_i = []; 
                for j=1:1:length(obj.subsystems)
                    A_i = [A_i, obj.subsystems(i).A{j}];
                    B_i = [B_i, obj.subsystems(i).B{j}];
                    C_i = [C_i, obj.subsystems(i).C{j}];
                    D_i = [D_i, obj.subsystems(i).D{j}];
                    E_i = [E_i, obj.subsystems(i).E{j}];
                    F_i = [F_i, obj.subsystems(i).F{j}];
                end
                A = [A; A_i]; B = [B; B_i]; C = [C; C_i]; D = [D; D_i]; E = [E; E_i]; F = [F; F_i]; 
                x = [x; obj.subsystems(i).x];
            end
            % storage incase
            obj.networkMatrices.A = A;
            obj.networkMatrices.B = B;
            obj.networkMatrices.C = C;
            obj.networkMatrices.D = D;
            obj.networkMatrices.E = E;
            obj.networkMatrices.F = F;
            obj.networkMatrices.x = x;
        end
        
        
        % update the network (mainly drawing stuff in the figure)
        function outputArg = updateNetwork(obj, figNum, timeVal)
            figure(figNum)
            
            % update targets
            thickness = 0.01; 
            scalingFactorx = 1; %1000
            scalingFactory = 1; %1000
            scalingFactoru = 1; %1000
            scalingFactorw = 1; %1000
            for i = 1:1:length(obj.subsystems)
                posX = obj.subsystems(i).location(1);
                posY = obj.subsystems(i).location(2);
                
                normx = norm(obj.subsystems(i).x); heightx = normx/scalingFactorx;
                normy = norm(obj.subsystems(i).y); heighty = normy/scalingFactory;
                normu = norm(obj.subsystems(i).u); heightu = normu/scalingFactoru;
                normw = norm(obj.subsystems(i).w); heightw = normw/scalingFactorw;
                
                if timeVal ~= 0
                    delete(obj.subsystems(i).graphicHandles(1));
                    delete(obj.subsystems(i).graphicHandles(2));
                    delete(obj.subsystems(i).graphicHandles(3));
                    delete(obj.subsystems(i).graphicHandles(4));
                end
                
                rx = rectangle('Position',[posX-3*thickness, posY, thickness, heightx],...
                    'FaceColor',[0.8 0 0 0.5],'EdgeColor','k','LineWidth',0.01);
%                 t = text(posX+0.02, posY+0.02,num2str(normx,4),'Color','b','FontSize',10);
                ry = rectangle('Position',[posX-2*thickness, posY, thickness, heighty],...
                    'FaceColor',[0 0.8 0 0.5],'EdgeColor','k','LineWidth',0.01);
                ru = rectangle('Position',[posX+thickness, posY, thickness, heightu],...
                    'FaceColor',[0 0 0.8 0.5],'EdgeColor','k','LineWidth',0.01);
                rw = rectangle('Position',[posX+2*thickness, posY, thickness, heightw],...
                    'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k','LineWidth',0.01);

                obj.subsystems(i).graphicHandles(1) = rx;
                obj.subsystems(i).graphicHandles(2) = ry;
                obj.subsystems(i).graphicHandles(3) = ru;
                obj.subsystems(i).graphicHandles(4) = rw; 
                
            end
        end
        
        
        % Ask subsystems to execute the computed update
        function outputArg = finishUpdate(obj)
            for i = 1:1:length(obj.subsystems)
                obj.subsystems(i).finishUpdate();
            end
        end
        
        
        
        function newIndexes = convertToNewPlacesList(~, originalIndexes, indexing)
            newIndexes = [];
            for i = 1:1:length(originalIndexes)
                iInd = originalIndexes(i); %original index
                newIndex = find(iInd==indexing); % new index according to the given indexing
                newIndexes = [newIndexes, newIndex];
            end
        end
        
        
        function minNewIndex = findMinNewPlace(~, originalIndexes, indexing, iTemp)
            minNewIndex = iTemp;
            for i = 1:1:length(indexing)
                if sum(indexing(i)==originalIndexes)==1
                    minNewIndex = i;
                    break
                end
            end
        end
        
        
        function out = loadLeastNeighbors(obj,indexing)
            % with respect to the given indexing scheme, re-assign the
            % least neighbors
            
            for i = 1:1:length(indexing)
                iInd = indexing(i); % this is an index of an original subsystem
                % obj.neighbors{iInd}      is the original indexes of the neighbors of iInd
                % obj.outNeighbors{iInd}   is the original indexes of the outneighbors of iInd
                % we need to map these sets of original indexes to the new
                % indexing and then pick the minimum new index.
                % this is what is being done in the following two lines
                obj.leastNeighbors(iInd) = obj.findMinNewPlace(obj.exclusiveNeighbors{iInd},indexing,i);       % 
                obj.leastOutNeighbors(iInd) = obj.findMinNewPlace(obj.exclusiveOutNeighbors{iInd},indexing,i); % original indexes of the neighbors of iInd
            end
            
        end
        
        
        % Get communication cost
        function comCost = getCommunicationCost(obj,indexing)
            
            obj.loadLeastNeighbors(indexing); % based on the new "indexing", find the new least neighbors (for \Gamma_i's and \Delta_i's)
            
            % Indexing is the list of preferred order of subsystems.
            comCost = 0;
            for i = 1:1:length(indexing)
                for j = 1:1:i-1
                    iInd = indexing(i);
                    jInd = indexing(j);
                    
                    % cost of having to communicate from j to i
                    alpha_ij = obj.distanceMatrix(jInd,iInd);
                    
                    % cost of having to send preliminary infomation to compute W_ij from j
                    % THis preliminary infomation is like terms A_ji, B_ji, ... 
                    beta_ij = 3*(obj.subsystems(jInd).dim_n)^2;
                    
                    % beta_ij needs to be sent to i from j only if i is a neighbor of j:
                    isiNeighborOfj = sum(iInd==obj.subsystems(jInd).neighbors)==1;
                    
                    % To send \tilde{W}_jj to i from j
                    gamma_ijj = (obj.subsystems(jInd).dim_n)^2;
                    
                    % \mathcal{L}_ij value
                    
                    Gamma_i = obj.leastNeighbors(iInd);
                    Delta_i = obj.leastOutNeighbors(iInd);
                    Gamma_j = obj.leastNeighbors(jInd);
                    Delta_j = obj.leastOutNeighbors(jInd);
                    L_ij = max(min(Gamma_i,Delta_i),min(Gamma_j,Delta_j));
                    
                    % cost of getting \tilde{W}_{jk} to i from j
                    sumgamma_ijk = 0;
                    if L_ij <= (j-1)
                        for k = L_ij:1:(j-1)
                            % To send \tilde{W}_{jk} to i from j
                            kInd = indexing(k);
                            gamma_ijk = obj.subsystems(jInd).dim_n*obj.subsystems(kInd).dim_n;
                            sumgamma_ijk = sumgamma_ijk + gamma_ijk;
                        end
                    end
                    
                    % Total cost at iteration i, to get thing relevent to iteration j
                    J_ij = alpha_ij*(beta_ij*isiNeighborOfj + gamma_ijj + sumgamma_ijk);
                    comCost = comCost + J_ij;
                    
                end
            end
        
        end
        
        
        function [bestIndexing, minCost, worstIndexing, maxCost] = findOptimumIndexing(obj)
            basicIndexing = 1:1:length(obj.subsystems);
            permutations = perms(basicIndexing);
            minCost = inf;
            maxCost = 0;
            L = size(permutations,1);
            for k = 1:1:L
%                 if rem(k,100)==0
%                     disp([num2str(k*100/L,3),"% is complete"])
%                 end
                
                cost = obj.getCommunicationCost(permutations(k,:));
                
                if cost < minCost
                    minCost = cost;
                    bestIndexing = permutations(k,:);
                end
                
                if cost > maxCost
                    maxCost = cost;
                    worstIndexing = permutations(k,:);
                end
            end
        end
        
        
        function outputArg = drawIndexing(obj,indexing)
            for i=1:1:length(indexing)
                iInd = indexing(i);
                obj.subsystems(iInd).drawIndex(i)
            end
        end
        
        
        
        % Assign given controller coefficients at subsystems
        function output = assignLocalControllers(obj,K)
            for i = 1:1:length(obj.subsystems)
                p = obj.subsystems(i).dim_p;
                for j = 1:1:length(obj.subsystems)
                    n = obj.subsystems(j).dim_n;
                    obj.subsystems(i).globalSFBLQRControllerGains{j} = K((i-1)*p+1:i*p,(j-1)*n+1:j*n);
                end
            end
        end
        
        
        % Assign a new A matrix at subsystes
        function changeNetworkParameters(obj,ANew)
            for i = 1:1:length(obj.subsystems)
                n_i = obj.subsystems(i).dim_n;
                for j = 1:1:length(obj.subsystems)
                    n_j = obj.subsystems(j).dim_n;
                    obj.subsystems(i).A{j} = ANew((i-1)*n_i+1:i*n_i,(j-1)*n_j+1:j*n_j);
                end
            end
        end
        
        
        function [Q,S,R] = getSomeQSRMatrices(obj,stringInput)
            dim_m = 0; % for Q
            dim_q = 0; % for R
            for i=1:1:length(obj.subsystems)
                dim_m = dim_m + obj.subsystems(i).dim_m;
                dim_q = dim_q + obj.subsystems(i).dim_q;
            end
            
            if isequal(stringInput,"strictly passive")
                nu = 1000000;
                rho = 1000000;               
                
                Q = -rho*eye(dim_m);
                S = 0.5*eye(dim_m,dim_q);
                R = -nu*eye(dim_q);
            elseif isequal(stringInput,"random")
                Q = 5*(rand(dim_m)-0.5); Q = Q + Q'; 
                S = 2*(rand(dim_m,dim_q)-0.5); 
                R = 5*(rand(dim_q)-0.5); R = R + R';
            else % passive
                Q = 0*eye(dim_m);
                S = 0.5*eye(dim_m,dim_q);
                R = 0*eye(dim_q);
            end
        end
        
        function output = storeQSRMatrices(obj,Q,S,R)
            for i = 1:1:length(obj.subsystems)
                m_i = obj.subsystems(i).dim_m;
                q_i = obj.subsystems(i).dim_q;
                for j = 1:1:length(obj.subsystems)
                    m_j = obj.subsystems(j).dim_m;
                    q_j = obj.subsystems(j).dim_q;
                    
                    obj.subsystems(i).dataToBeDistributed.Q{j} = Q((i-1)*m_i+1:i*m_i,(j-1)*m_j+1:j*m_j);
                    obj.subsystems(i).dataToBeDistributed.S{j} = S((i-1)*m_i+1:i*m_i,(j-1)*q_j+1:j*q_j);
                    obj.subsystems(i).dataToBeDistributed.R{j} = R((i-1)*q_i+1:i*q_i,(j-1)*q_j+1:j*q_j);
                    
                end
            end
            % Store at network level as well
            obj.networkMatrices.Q = Q;
            obj.networkMatrices.S = S;
            obj.networkMatrices.R = R;
        end 
        
        
        
        %% Centralized stability and dissipativity tests
        
        % Centralized stability test
        function output = centralizedStabilityTest(obj)
            A = obj.networkMatrices.A;
            setlmis([]);  % To initialize the LMI description
            P = lmivar(1,[size(A,1), 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 
            lmiterm([-1, 1, 1, P],-1,A,'s');
            lmiterm([-2, 1, 1, P],1,1); % defines -P<0
            lmisys = getlmis;
            [tmin,~] = feasp(lmisys); % Solve the LMI system
            output = tmin <= 0; % feasible if this is satisfied
        end
        
            
        % Centralized dissipativity test
        function output = centralizedQSRDissipativityTest(obj)
            A = obj.networkMatrices.A;
            E = obj.networkMatrices.E; % instead of B
            C = obj.networkMatrices.C; 
            F = obj.networkMatrices.F; % instead of D
            Q = obj.networkMatrices.Q;
            S = obj.networkMatrices.S;
            R = obj.networkMatrices.R;
            
            setlmis([]);  % To initialize the LMI description
            P = lmivar(1,[size(A,1), 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 

            % W = [W_ii1, W_2; W_3, W_4]
            % W_1 = -A^T P - P A + C^\T Q C 
            lmiterm([-1, 1, 1, P],-1,A,'s');
            W1 = C'*Q'*C;
            lmiterm([-1, 1, 1, 0],W1);
            % W_2 = -PE + C^\T S + C^\T Q F
            lmiterm([-1, 1, 2, P],-1,E);
            W2 = C'*S + C'*Q*F;
            lmiterm([-1, 1, 2, 0],W2);
            % W_3 = (-PE + C^\T S + C^\T Q F)^\T 
            lmiterm([-1, 2, 1, P],-E',1);
            lmiterm([-1, 2, 1, 0],W2');
            % W_4 = F^\T Q F + (F^\T S + S^\T F) + R
            W4 = F'*Q*F + (F'*S + S'*F) + R;
            lmiterm([-1, 2, 2, 0],W4');
            
            % P>0
            lmiterm([-2, 1, 1, P],1,1); % defines -P<0
            lmisys = getlmis;
            [tmin,~] = feasp(lmisys); % Solve the LMI system
            output = tmin <= 0; % strictly feasible if this is satisfied
        end
        
        
        
        %% Decentralized stability and dissipativity tests
        
        
        % Check stability of the networked system
        function output = checkStability(obj,indexing,typeVal)
            
            if isempty(indexing)
                indexing = [1:1:length(obj.subsystems)];
            end
            
            for i = 1:1:length(indexing)
                iInd = indexing(i);
                previousSubsystems = indexing(1:i-1);
                if typeVal==1
                    isStable = obj.subsystems(iInd).checkStability(previousSubsystems,obj.subsystems);
                else
                    isStable = obj.subsystems(iInd).checkStability2(previousSubsystems,obj.subsystems);
                end
                if ~isStable
                    output = false;
                    break
                else
                    output = true;
                end
            end
            
        end
        
        
        % Check QSR-dissipativity of the networked system
        function output = checkQSRDissipativity(obj,indexing)
            
            if isempty(indexing)
                indexing = [1:1:length(obj.subsystems)];
            end
            
            for i = 1:1:length(indexing)
                iInd = indexing(i);
                previousSubsystems = indexing(1:i-1);                
                isQSRDissipative = obj.subsystems(iInd).checkQSRDissipativity(previousSubsystems, obj.subsystems);
                if ~isQSRDissipative
                    output = false;
                    break
                else
                    output = true;
                end
            end
            
        end
        
        
        
        
        
        %% Centralized control design
        
        % Stabilizing state-feedback controller design
        function [K, isStabilizable] = designGlobalStabilizingSFBControllers(obj)
            A = obj.networkMatrices.A;
            B = obj.networkMatrices.B;
            
            setlmis([]);  % To initialize the LMI description
            Q = lmivar(1,[size(A,1), 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 
            L = lmivar(2,[size(B,2),size(A,1)]);
            
            % We need Q>0 and -AQ-BL-QA'-L'B'>0 where P = Q^{-1} and K = LQ^{-1}
            lmiterm([-1, 1, 1, Q],-1,A','s'); % defines -AQ-QA'>0
            lmiterm([-1, 1, 1, L],-B,1,'s');   % defines -BL-L'B'>0
            lmiterm([-2, 1, 1, Q],1,1);       % defines -Q<0
            lmisys = getlmis;
            [tmin,sol] = feasp(lmisys); % Solve the LMI system
            isStabilizable = tmin <= 0; % feasible if this is satisfied
            Q = dec2mat(lmisys,sol,Q); 
            L = dec2mat(lmisys,sol,L); 
            K = L/Q;
        end
        
        % Stabilizing state-feedback controller design
        function [K, isDissipative] = designGlobalDissipatingSFBControllers(obj)
            A = obj.networkMatrices.A;
            B = obj.networkMatrices.B;
            E = obj.networkMatrices.E; % instead of B
            C = obj.networkMatrices.C; 
            F = obj.networkMatrices.F; % instead of D
            Q = obj.networkMatrices.Q;
            S = obj.networkMatrices.S;
            R = obj.networkMatrices.R;
            
            setlmis([]);  % To initialize the LMI description
            L = lmivar(1,[size(A,1), 1]); % L = P^{-1} > 0 is a variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 
            M = lmivar(2,[size(B,2), size(A,1)]); % M = KP^{-1} is a free variable
            
            % W = [W_11, W_12, W_13; W_21, W_22, W_23; W_31, W_32, W_34]
            
            % W_11 =  -LA^T -AL - M^T B^T - BM
            lmiterm([-1, 1, 1, L],-1,A','s');
            lmiterm([-1, 1, 1, M],-B,1,'s');
            % W_12 = -E + L C^T S
            lmiterm([-1, 1, 2, L],1,C'*S);
            lmiterm([-1, 1, 2, 0],-E);
            % W_13 = LC^T
            lmiterm([-1, 1, 3, L],1,C');
            
            % W_21 = -E^T + S^T C L
            lmiterm([-1, 2, 1, L],S'*C,1);
            lmiterm([-1, 2, 1, 0],-E');
            % W_22 = F^T S + S^T F + R
            lmiterm([-1, 2, 2, 0],F'*S+S'*F+R);
            % W_23 = F^T
            lmiterm([-1, 2, 3, 0],F');
            
            % W_31 = CL
            lmiterm([-1, 3, 1, L],C,1);
            % W_32 = F
            lmiterm([-1, 3, 2, 0],F);
            % W_33 = Q^{-1}
            lmiterm([-1, 3, 3, 0],inv(Q));
            
            % L>0
            lmiterm([-2, 1, 1, L],1,1); % defines -P<0
            lmisys = getlmis;
            [tmin,sol] = feasp(lmisys); % Solve the LMI system
            isDissipative = tmin <= 0; % strictly feasible if this is satisfied
            L = dec2mat(lmisys,sol,L); 
            M = dec2mat(lmisys,sol,M); 
            K = M/L;
        end
        
        
        
        %% Decentralized control design
        
        % Stabilizing state-feedback controller design
        function [K, isStabilizable] = designLocalStabilizingSFBControllers(obj,indexing)
        
            if isempty(indexing)
                indexing = [1:1:length(obj.subsystems)];
            end
            
            for i = 1:1:length(indexing)
                iInd = indexing(i);
                previousSubsystems = indexing(1:i-1);                
                [isStabilizable,K_ii,K_ijVals,K_jiVals] = obj.subsystems(iInd).designLocalStabilizingSFBControllers(previousSubsystems, obj.subsystems);
                
                K{iInd,iInd} = K_ii;
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    obj.subsystems(jInd).localStabilizingSFBControllerGains{iInd} = K_jiVals{jInd};
                    
                    K{iInd,jInd} = K_ijVals{jInd};
                    K{jInd,iInd} = K_jiVals{jInd};
                end
                
                if ~isStabilizable
                    break
                end
                
            end
            
            % Collect all the coefficients
        end
        
        
        % Stabilizing state-feedback controller design
        function [K, isDissipative] = designLocalDissipatingSFBControllers(obj,indexing)
        
            if isempty(indexing)
                indexing = [1:1:length(obj.subsystems)];
            end
            
            K = [];
            for i = 1:1:length(indexing)
                iInd = indexing(i);
                previousSubsystems = indexing(1:i-1);                
                [isDissipative,K_ii,K_ijVals,K_jiVals] = obj.subsystems(iInd).designLocalDissipatingSFBControllers(previousSubsystems, obj.subsystems);
                
                K{iInd,iInd} = K_ii;
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    obj.subsystems(jInd).localStabilizingSFBControllerGains{iInd} = K_jiVals{jInd};
                    
                    K{iInd,jInd} = K_ijVals{jInd};
                    K{jInd,iInd} = K_jiVals{jInd};
                end
                
                if ~isDissipative
                    break
                end
            end
                    
        end
        
        
       
    end
end

