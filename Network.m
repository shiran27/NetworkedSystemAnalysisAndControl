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
        function obj = loadARandomNetwork(obj,numOfSubsystems,dimentionOfSpace,sizeOfSpace,communicationRadius,dims,diags)
            
            % Generate subsystems' list
            subsystemLocations = sizeOfSpace*rand(numOfSubsystems,dimentionOfSpace);
            
            % Initializing subsystems: 
            % Subsystem-spaecific parameters will be loaded inside
            obj.subsystems = [];
            for i = 1:1:numOfSubsystems
                newSubsystem = Subsystem(i, subsystemLocations(i,:),dims{i});
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
                % diags = nature of [A,B,C,D,E,F,G,H,J], 1 diag, 0 general, -1 zero 
%                 obj.subsystems(i).loadParameters(obj.subsystems);         % Completely random
                obj.subsystems(i).loadStableParameters(obj.subsystems,diags);   % Lets use this in the journal paper
%                 obj.subsystems(i).loadPassiveParameters(obj.subsystems);  % this was used in the Med conf paper
                
%                 obj.subsystems(i).designLocalSFBLQRControllerGains(); % find local state feedback LQR controller gains
                
            end
            
            
        end
        
        % function to load a custom networked system
        function obj = loadTheCustomNetwork(obj)
            
            % Subsystem locations (for plotting purposes only) Nx2 matrix
%             subsystemLocations = [0.5,0.5; 0.2, 0.2; 0.8, 0.2; 0.8,0.8; 0.2,0.8];
            subsystemLocations = [0.25,0.5; 0.5, 0.5; 0.75, 0.5; 0.375,0.25];
            numOfSubsystems = size(subsystemLocations,1);
            
            % Initializing subsystems: 
            % Subsystem-spaecific parameters will be loaded inside
            obj.subsystems = [];
            for i = 1:1:numOfSubsystems
                dims.n = 4; % x
                dims.p = 2; % u
                dims.q = 2; % w  
                dims.m = 2; % y
                newSubsystem = Subsystem(i, subsystemLocations(i,:),dims);
                obj.subsystems = [obj.subsystems, newSubsystem];
            end
            
            % Inter-subsystem connections: defined by a list of edges Nx2
%             edgeList = [1,2;2,1;1,3;3,1;1,4;4,1;1,5;5,1;2,3;3,4;4,5;5,2];
            edgeList = [1,2;2,1;2,3;3,2;1,4;4,1;4,2;2,4];
        
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
%                 obj.subsystems(i).designLocalSFBLQRControllerGains(); % find local state feedback LQR controller gains     
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
        function outputArg = drawNetwork(obj, figNum, args)
            
            sizeOfSpace = 1;
            legendOn = args(1);
            if length(args)>1
                axisLimits = args(2:5)
            else
                axisLimits = [0,sizeOfSpace,0,sizeOfSpace];
            end
            
            figure(figNum)
            
            grid on
            axis equal
            axis(axisLimits);
            axis off
            
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
            
            if legendOn
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
        end
        
        
        function output = shiftLocations(obj)
            centerX = 0;
            centerY = 0;
            for i = 1:1:length(obj.subsystems)
                centerX = centerX + obj.subsystems(i).location(1);
                centerY = centerY + obj.subsystems(i).location(2);
            end
            
            xShift = 0.5 - centerX/length(obj.subsystems);
            yShift = 0.5 - centerY/length(obj.subsystems);
            
            for i = 1:1:length(obj.subsystems)
                obj.subsystems(i).location(1) = obj.subsystems(i).location(1) + xShift;
                obj.subsystems(i).location(2) = obj.subsystems(i).location(2) + yShift;
            end
            
            for i = 1:1:length(obj.edges)
                obj.edges(i).locations(:,1) = obj.edges(i).locations(:,1) + xShift;
                obj.edges(i).locations(:,2) = obj.edges(i).locations(:,2) + yShift;
            end
        end
        
        % Generate the complete system.
        function [A,B,C,D,E,F,G,H,J,x] = getNetworkMatrices(obj)
            A = []; B = []; C = []; D = []; E = []; F = [];  G = [];  H = [];  J = []; x = [];
            for i=1:1:length(obj.subsystems)
                A_i = []; B_i = []; C_i = []; D_i = []; E_i = []; F_i = []; G_i = []; H_i = []; J_i = [];
                for j=1:1:length(obj.subsystems)
                    A_i = [A_i, obj.subsystems(i).A{j}];
                    B_i = [B_i, obj.subsystems(i).B{j}];
                    C_i = [C_i, obj.subsystems(i).C{j}];
                    D_i = [D_i, obj.subsystems(i).D{j}];
                    E_i = [E_i, obj.subsystems(i).E{j}];
                    F_i = [F_i, obj.subsystems(i).F{j}];
                    G_i = [G_i, obj.subsystems(i).G{j}];
                    H_i = [H_i, obj.subsystems(i).H{j}];
                    J_i = [J_i, obj.subsystems(i).J{j}];
                end
                A = [A; A_i]; B = [B; B_i]; C = [C; C_i]; D = [D; D_i]; E = [E; E_i]; F = [F; F_i];  G = [G; G_i]; H = [H; H_i]; J = [J; J_i];
                x = [x; obj.subsystems(i).x];
            end
            % storage incase
            obj.networkMatrices.A = A;
            obj.networkMatrices.B = B;
            obj.networkMatrices.C = C;
            obj.networkMatrices.D = D;
            obj.networkMatrices.E = E;
            obj.networkMatrices.F = F;
            obj.networkMatrices.G = G;
            obj.networkMatrices.H = H;
            obj.networkMatrices.J = J;
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
                    alpha_ij = obj.distanceMatrix(jInd,iInd)>1.5;
                    
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
        
        
        function [bestIndexing, minCost, worstIndexing, maxCost, basicIndexingCost] = findOptimumIndexing(obj)
            basicIndexing = 1:1:length(obj.subsystems);
            basicIndexingCost = obj.getCommunicationCost(basicIndexing);
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
        
        
        function outputArg = drawIndexing(obj,indexing,angleVal,colorVal)
            for i=1:1:length(indexing)
                iInd = indexing(i);
                obj.subsystems(iInd).drawIndex(i,angleVal,colorVal)
            end
        end
        
        
        
        % Assign given controller coefficients at subsystems
        function output = assignLocalControllers(obj,K,typeVal)
            for i = 1:1:length(obj.subsystems)
                p = obj.subsystems(i).dim_p;
                for j = 1:1:length(obj.subsystems)
                    n = obj.subsystems(j).dim_n;
                    if any(obj.subsystems(i).neighbors==j)
                        Kblock = K((i-1)*p+1:i*p,(j-1)*n+1:j*n);
                    else
                        Kblock = zeros(p,n);
                    end
                    
                    if isequal(typeVal,'globalSFBLQR')
                        obj.subsystems(i).globalSFBLQRControllerGains{j} = Kblock;
                    elseif isequal(typeVal,'globalSFBLMI')
                        obj.subsystems(i).globalSFBLMIControllerGains{j} = Kblock;
                    end
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
        
        
        function [Q,S,R] = getSomeQSRMatrices(obj,dissFrom,dissTo,dissType,dissArgs)
            
%             dissFrom = 'w';    % Options: 'u', 'w'
%             dissTo = 'y';      % Options: 'y', 'z',
%             dissType = 'strictly passive';     % Options: strictly passive([IFP, OFP]), passive, L2G(gamma), random 
%             dissArgs = [0.001,1]
%             [Q,S,R] = network.getSomeQSRMatrices(dissFrom,dissTo,dissType,dissArgs) 
            
            dim_in = 0; % For R matrix
            dim_out = 0; % For Q matrix
            for i=1:1:length(obj.subsystems)
                
                if isequal(dissFrom,'u')
                    dim_in = dim_in + obj.subsystems(i).dim_p;
                elseif  isequal(dissFrom,'w')
                    dim_in = dim_in + obj.subsystems(i).dim_q;
                end
                
                if isequal(dissTo,'y')
                    dim_out = dim_out + obj.subsystems(i).dim_m;
                elseif isequal(dissTo,'z')
                    dim_out = dim_out + obj.subsystems(i).dim_l;
                end 
                
            end
            
            
            if isequal(dissType,'strictly passive')
                nu = dissArgs(1);   % IFP
                rho = dissArgs(2);  % OFP             
                Q = -rho*eye(dim_out);
                S = 0.5*eye(dim_out,dim_in);
                R = -nu*eye(dim_in);
            elseif isequal(dissType,'L2G')
                gamma = dissArgs(1); % L2G
                Q = -(1/gamma)*eye(dim_out);
                S = zeros(dim_out,dim_in);
                R = -gamma*eye(dim_in);
            elseif isequal(dissType,'passive')
                Q = zeros(dim_out);
                S = 0.5*eye(dim_out,dim_in);
                R = zeros(dim_in);
            elseif isequal(dissType,'random')
                Q = 5*(rand(dim_out)-0.5); Q = Q*Q.'; % Q = Q + Q'; 
                S = 2*(rand(dim_out,dim_in)-0.5); 
                R = 5*(rand(dim_in)-0.5); R = R*R.'; % R = R + R';
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
        
        
        
        %% Centralized stability based results
        
        % Centralized stability analysis
        function isStable = centralizedStabilityAnalysis(obj,solverOptions)
            A = obj.networkMatrices.A;
            % n = size(A,1);
            
            % P = sdpvar(n,n)
            % P = diag(sdpvar(n,1)); % Stability analysis with diagonal matrix
            P = [];
            for i = 1:1:length(obj.subsystems)
                n_i = obj.subsystems(i).dim_n;
                P = blkdiag(P,sdpvar(n_i,n_i));
            end
            
            con1 = P >= 0;
            con2 = A'*P + P*A <= 0;
            sol = optimize([con1,con2],[],solverOptions);
            isStable = sol.problem==0;
%             eigs = eig(A)
%             PVal = value(P)
        end
        
        
        % Centralized FSF stabilization
        function [K, isFSFStabilizable] = centralizedFSFStabilization(obj,solverOptions)
            A = obj.networkMatrices.A;
            B = obj.networkMatrices.B;
            % n = size(A,1);
            % p = size(B,2);
            
            % M = sdpvar(n,n)
            % M = diag(sdpvar(n,1));
            % L = sdpvar(p,n);
            M = [];
            L = [];
            for i = 1:1:length(obj.subsystems)
                n_i = obj.subsystems(i).dim_n;
                p_i = obj.subsystems(i).dim_p;
                M = blkdiag(M,sdpvar(n_i,n_i));
                L_i = [];
                for j = 1:1:length(obj.subsystems)
                    n_j = obj.subsystems(j).dim_n;
                    if any(obj.neighbors{i}==j)
                        L_i = [L_i,sdpvar(p_i,n_j)];
                    else
                        L_i = [L_i,zeros(p_i,n_j)];
                    end
                end
                L = [L;L_i];
            end
            
            con1 = M >= 0;
            con2 = M*A' + A*M + L'*B' + B*L <= 0;
            sol = optimize([con1,con2],[],solverOptions);
            isFSFStabilizable = sol.problem==0;
            MVal = value(M);
            LVal = value(L);
            K = LVal/MVal;
%             eigs = eig(A+B*K)
            
        end
        
        
        % Centralized stable observer design
        function [L, isObserverStable] = centralizedStableObserverDesign(obj,solverOptions)
            A = obj.networkMatrices.A;
            C = obj.networkMatrices.C;
            % n = size(A,1);
            % m = size(C,1);
            
            % P = sdpvar(n,n)
            % P = diag(sdpvar(n,1));
            % K = sdpvar(n,m);
            
            P = [];
            K = [];
            for i = 1:1:length(obj.subsystems)
                n_i = obj.subsystems(i).dim_n;
                P = blkdiag(P,sdpvar(n_i,n_i));
                K_i = [];
                for j = 1:1:length(obj.subsystems)
                    m_j = obj.subsystems(j).dim_m;
                    if any(obj.neighbors{i}==j)
                        K_i = [K_i,sdpvar(n_i,m_j)];
                    else
                        K_i = [K_i,zeros(n_i,m_j)];
                    end
                end
                K = [K;K_i];
            end
            
            con1 = P >= 0;
            con2 = A'*P + P*A - C'*K' - K*C <= 0;
            sol = optimize([con1,con2],[],solverOptions);
            isObserverStable = sol.problem==0;
            PVal = value(P);
            KVal = value(K);
            L = PVal\KVal;
%             eigs = eig(A-L*C)
        end
        
        
        % Centralized DOF stabilization
        function [Ac,Bc,Cc,Dc, isDOFStabilizable] = centralizedDOFStabilization(obj,solverOptions)
            A = obj.networkMatrices.A;
            B = obj.networkMatrices.B;
            C = obj.networkMatrices.C;
            D = obj.networkMatrices.D;
            if sum(sum(abs(D)))~=0
                disp('System matrix D is not null.')
                return
            end
            n = size(A,1);
            % p = size(B,2);
            % q = size(C,1);
            
            % Y = diag(sdpvar(n,1));
            % X = diag(sdpvar(n,1));
            % An = sdpvar(n,n,'full');
            % Bn = sdpvar(n,q);
            % Cn = sdpvar(p,n);
            % Dn = sdpvar(p,q);
            
            X = []; Y = [];
            An = []; Bn = []; Cn = []; Dn = [];
            for i = 1:1:length(obj.subsystems)
                n_i = obj.subsystems(i).dim_n;
                p_i = obj.subsystems(i).dim_p;
                X = blkdiag(X,sdpvar(n_i,n_i));
                Y = blkdiag(Y,sdpvar(n_i,n_i));
                An_i = []; Bn_i = []; Cn_i = []; Dn_i = [];
                for j = 1:1:length(obj.subsystems)
                    n_j = obj.subsystems(j).dim_n;
                    p_j = obj.subsystems(j).dim_p;
                    m_j = obj.subsystems(j).dim_m;                    
                    if any(obj.neighbors{i}==j)
                        An_i = [An_i, sdpvar(n_i,n_j)];
                        Bn_i = [Bn_i, sdpvar(n_i,m_j)];
                        Cn_i = [Cn_i, sdpvar(p_i,n_j)];
                        Dn_i = [Dn_i, sdpvar(p_i,m_j)];
                    else
                        An_i = [An_i, zeros(n_i,n_j)];
                        Bn_i = [Bn_i, zeros(n_i,m_j)];
                        Cn_i = [Cn_i, zeros(p_i,n_j)];
                        Dn_i = [Dn_i, zeros(p_i,m_j)];
                    end
                end
                An = [An; An_i]; Bn = [Bn; Bn_i]; Cn = [Cn; Cn_i]; Dn = [Dn; Dn_i];
            end
            
            I = eye(n);
            con1 = X >= 0;
            con2 = Y >= 0;
            con3 = [Y, I; I, X] >= 0;
            c1 = A*Y+B*Cn;
            c2 = A+B*Dn*C+An';
            c3 = X*A+Bn*C;
            con4 = [-c1-c1', -c2; -c2', -c3-c3']>=0;
            
            sol = optimize([con1,con2,con3,con4],[],solverOptions);
            isDOFStabilizable = sol.problem==0;
            
            XVal = value(X);
            YVal = value(Y);
            AnVal = value(An);
            BnVal = value(Bn);
            CnVal = value(Cn);
            DnVal = value(Dn);
            [M,N] = lu(I-XVal*YVal); % controller parameters
            N = N';
            Ac = M\(AnVal-BnVal*C*YVal-XVal*B*CnVal-XVal*(A-B*DnVal*C)*YVal)/(N'); 
            Bc = M\(BnVal-XVal*B*DnVal);
            Cc = (CnVal-DnVal*C*YVal)/(N');
            Dc = DnVal;
            Abar = [A+B*Dc*C, B*Cc; Bc*C, Ac]; %closed loop system
%             eigs = eig(Abar)
        end
        
 
           
        %% Centralized dissipativity based results
        
        % Centralized dissipativity analysis
        function isDissipative = centralizedDissipativityAnalysis(obj,dissFrom,dissTo,Q,S,R,solverOptions)
            A = obj.networkMatrices.A;
            if isequal(dissFrom,'u')  
                B = obj.networkMatrices.B;
                if isequal(dissTo,'y')
                    C = obj.networkMatrices.C;
                    D = obj.networkMatrices.D;
                elseif isequal(dissTo,'z')
                    C = obj.networkMatrices.G;
                    D = obj.networkMatrices.H;
                end                
            elseif isequal(dissFrom,'w')
                B = obj.networkMatrices.E;
                if isequal(dissTo,'y')
                    C = obj.networkMatrices.C;
                    D = obj.networkMatrices.F;                    
                elseif isequal(dissTo,'z')
                    C = obj.networkMatrices.G;
                    D = obj.networkMatrices.J;
                end
            end
            % n = size(A,1);
            
            % P = sdpvar(n,n)
            % P = diag(sdpvar(n,1)); % Stability analysis with diagonal matrix
            P = [];
            for i = 1:1:length(obj.subsystems)
                n_i = obj.subsystems(i).dim_n;
                P = blkdiag(P,sdpvar(n_i,n_i));
            end
            
            con1 = P >= 0;
            c1 = -A'*P-P*A;
            c2 = -P*B+C'*S;
            c3 = D'*S+S'*D+R;
            if all(Q(:)==0)
                con2 = [c1, c2; c2', c3] >= 0;
            else
                con2 = [c1, c2, C'; c2', c3, D'; C, D, -inv(Q)] >= 0;
            end
            sol = optimize([con1,con2],[],solverOptions);
            isDissipative = sol.problem==0;
%             eigs = eig(A)
        end
        
        
        % Centralized FSF dissipativation
        function [K, isFSFDissipative] = centralizedFSFDissipativation(obj,dissFrom,dissTo,Q,S,R,solverOptions)
            A = obj.networkMatrices.A;
            B = obj.networkMatrices.B;
            E = obj.networkMatrices.E;
            if isequal(dissFrom,'w') % the only possibility under FSF
                if isequal(dissTo,'y')
                    C = obj.networkMatrices.C;
                    D = obj.networkMatrices.D;
                    F = obj.networkMatrices.F;                    
                elseif isequal(dissTo,'z')
                    C = obj.networkMatrices.G;
                    D = obj.networkMatrices.H;
                    F = obj.networkMatrices.J;
                end
            end
            if sum(sum(abs(D)))~=0
                disp('System matrix D is not null!')
                return
            end
            % n = size(A,1);
            
            % M = sdpvar(n,n)
            % M = diag(sdpvar(n,1)); % Stability analysis with diagonal matrix
            % L = sdpvar(1,n);
            M = [];
            L = [];
            for i = 1:1:length(obj.subsystems)
                n_i = obj.subsystems(i).dim_n;
                p_i = obj.subsystems(i).dim_p;
                M = blkdiag(M,sdpvar(n_i,n_i));
                L_i = [];
                for j = 1:1:length(obj.subsystems)
                    n_j = obj.subsystems(j).dim_n;
                    if any(obj.neighbors{i}==j)
                        L_i = [L_i,sdpvar(p_i,n_j,'full')];
                    else
                        L_i = [L_i,zeros(p_i,n_j)];
                    end
                end
                L = [L;L_i];
            end
            
            con1 = M >= 0;
            c1 = A*M + B*L;
            c2 = -E+M*C'*S;
            c3 = F'*S+S'*F+R;
            if all(Q(:)==0)
                con2 = [-c1-c1', c2; c2', c3] >= 0;
            else
                con2 = [-c1-c1', c2, M*C'; c2', c3, F'; C*M', F, -inv(Q)] >= 0;
            end
            sol = optimize([con1,con2],[],solverOptions);
            isFSFDissipative = sol.problem==0;
            MVal = value(M);
            LVal = value(L);
            K = LVal/MVal;
%             eigs = eig(A+B*K)
        end
        
        
        % Centralized dissipative observer design
        function [L, isObserverDissipative] = centralizedDissipativeObserverDesign(obj,dissFrom,dissTo,Q,S,R,solverOptions)
            A = obj.networkMatrices.A;
            C = obj.networkMatrices.C;
            E = obj.networkMatrices.E;
            F = obj.networkMatrices.F;
            if isequal(dissFrom,'w') % the only possibility under FSF
                if isequal(dissTo,'z')
                    G = obj.networkMatrices.G;
                    J = obj.networkMatrices.J;
                end
            end
            % n = size(A,1);
            % m = size(C,1);
            
            
            % P = sdpvar(n,n)
            % P = diag(sdpvar(n,1));
            % K = sdpvar(n,m);
            P = [];
            K = [];
            for i = 1:1:length(obj.subsystems)
                n_i = obj.subsystems(i).dim_n;
                P = blkdiag(P,sdpvar(n_i,n_i));
                K_i = [];
                for j = 1:1:length(obj.subsystems)
                    m_j = obj.subsystems(j).dim_m;
                    if any(obj.neighbors{i}==j)
                        K_i = [K_i,sdpvar(n_i,m_j,'full')];
                    else
                        K_i = [K_i,zeros(n_i,m_j)];
                    end
                end
                K = [K;K_i];
            end
            
            
            con1 = P >= 0;
            c1 = P*A-K*C;
            c2 = -P*E+K*F+G'*S;
            c3 = J'*S+S'*J+R;
            if all(Q(:)==0)
                con2 = [-c1-c1', c2; c2', c3] >= 0;
            else
                con2 = [-c1-c1', c2, G'; c2', c3, J'; G, J, -inv(Q)] >= 0;
            end            
            sol = optimize([con1,con2],[],solverOptions);
            isObserverDissipative = sol.problem==0;
            PVal = value(P);
            KVal = value(K);
            L = PVal\KVal;
%             eigs = eig(A-L*C)           
            
        end
        
        
        % Centralized DOF dissipativation
        function [Ac,Bc,Cc,Dc, isDOFDissipative] = centralizedDOFDissipativation(obj,dissFrom,dissTo,Q,S,R,solverOptions)
            A = obj.networkMatrices.A;
            B = obj.networkMatrices.B;
            C = obj.networkMatrices.C;
            D = obj.networkMatrices.D;
            E = obj.networkMatrices.E;
            F = obj.networkMatrices.F;
            if isequal(dissFrom,'w') % the only possibility under DOF
                if isequal(dissTo,'y')
                    G = obj.networkMatrices.C;
                    H = obj.networkMatrices.D;
                    J = obj.networkMatrices.F;                    
                elseif isequal(dissTo,'z')
                    G = obj.networkMatrices.G;
                    H = obj.networkMatrices.H;
                    J = obj.networkMatrices.J;
                end
            end
            if sum(sum(abs(D)))~=0
                disp('System matrix D is not null!')
                return
            end
            n = size(A,1);
            % p = size(B,2);
            % q = size(C,1);
            
            % Y = diag(sdpvar(n,1));
            % X = diag(sdpvar(n,1));
            % An = sdpvar(n,n,'full');
            % Bn = sdpvar(n,q);
            % Cn = sdpvar(p,n);
            % Dn = sdpvar(p,q);            
            X = []; Y = [];
            An = []; Bn = []; Cn = []; Dn = [];
            for i = 1:1:length(obj.subsystems)
                n_i = obj.subsystems(i).dim_n;
                p_i = obj.subsystems(i).dim_p;
                X = blkdiag(X,sdpvar(n_i,n_i));
                Y = blkdiag(Y,sdpvar(n_i,n_i));
                An_i = []; Bn_i = []; Cn_i = []; Dn_i = [];
                for j = 1:1:length(obj.subsystems)
                    n_j = obj.subsystems(j).dim_n;
                    p_j = obj.subsystems(j).dim_p;
                    m_j = obj.subsystems(j).dim_m;                    
                    if any(obj.neighbors{i}==j)
                        An_i = [An_i, sdpvar(n_i,n_j,'full')];
                        Bn_i = [Bn_i, sdpvar(n_i,m_j,'full')];
                        Cn_i = [Cn_i, sdpvar(p_i,n_j,'full')];
                        Dn_i = [Dn_i, sdpvar(p_i,m_j,'full')];
                    else
                        An_i = [An_i, zeros(n_i,n_j)];
                        Bn_i = [Bn_i, zeros(n_i,m_j)];
                        Cn_i = [Cn_i, zeros(p_i,n_j)];
                        Dn_i = [Dn_i, zeros(p_i,m_j)];
                    end
                end
                An = [An; An_i]; Bn = [Bn; Bn_i]; Cn = [Cn; Cn_i]; Dn = [Dn; Dn_i];
            end
            
            
            I = eye(n);
            con1 = X >= 0;
            con2 = Y >= 0;
            con3 = [Y, I; I, X] >= 0;
            c1 = -A*Y-B*Cn;
            c2 = -A-B*Dn*C-An';
            c3 = -E-B*Dn*F+(Y*G'+Cn'*H')*S;
            c4 = Y*G'+Cn'*H';
            c5 = -X*A-Bn*C;
            c6 = -X*E-Bn*F+(G'+C'*Dn'*H')*S;
            c7 = G'+C'*Dn'*H';
            c8 = (J'+F'*Dn'*H')*S+S'*(J+H*Dn*F)+R;
            c9 = J'+F'*Dn'*H';
            if all(Q(:)==0)
                con4 = [c1+c1', c2, c3; c2', c5+c5', c6; c3', c6', c8]>=0;
            else
                con4 = [c1+c1', c2, c3, c4; c2', c5+c5', c6, c7; c3', c6', c8, c9; c4', c7', c9', -inv(Q)]>=0;
            end 
            sol = optimize([con1,con2,con3,con4],[],solverOptions);
            isDOFDissipative = sol.problem==0;
            XVal = value(X);
            YVal = value(Y);
            AnVal = value(An);
            BnVal = value(Bn);
            CnVal = value(Cn);
            DnVal = value(Dn);
            [M,N] = lu(I-XVal*YVal); % controller parameters
            N = N';
            Ac = M\(AnVal-BnVal*C*YVal-XVal*B*CnVal-XVal*(A-B*DnVal*C)*YVal)/(N');
            Bc = M\(BnVal-XVal*B*DnVal);
            Cc = (CnVal-DnVal*C*YVal)/(N');
            Dc = DnVal;
            Abar = [A+B*Dc*C, B*Cc; Bc*C, Ac]; %closed loop system
%             Bbar = [E+B*Dc*F;Bc*F];
%             Cbar = [G+H*Dc*C, H*Cc];
%             Dbar = [J+H*Dc*F];
%             eigs = eig(Abar)
        end
        
        
              
        
        
        %% Decentralized stability related concepts
        
        % Decentralized stability analysis
        function output = decentralizedStabilityAnalysis(obj,indexing,solverOptions)
            
            if isempty(indexing)
                indexing = [1:1:length(obj.subsystems)];
            end
            
            for i = 1:1:length(indexing)
                iInd = indexing(i);
                previousSubsystems = indexing(1:i-1);                
                isStable = obj.subsystems(iInd).stabilityAnalysis(previousSubsystems,obj.subsystems,solverOptions);
                
                if ~isStable
                    output = false;
                    break
                else
                    output = true;
                end
            end
            
        end
        
        
        % Decentralized FSF stabilization
        function [K, isStabilizable] = decentralizedFSFStabilization(obj,indexing,solverOptions)
        
            if isempty(indexing)
                indexing = [1:1:length(obj.subsystems)];
            end
            
            for i = 1:1:length(indexing)
                iInd = indexing(i);
                previousSubsystems = indexing(1:i-1);                
                [isStabilizable,K_ii,K_ijVals,K_jiVals] = obj.subsystems(iInd).FSFStabilization(previousSubsystems, obj.subsystems, solverOptions);
                
                K{iInd,iInd} = K_ii;
                obj.subsystems(iInd).controllerGains.decenFSFStabCont{iInd} = K_ii;
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    obj.subsystems(jInd).controllerGains.decenFSFStabCont{iInd} = K_jiVals{jInd};
                    obj.subsystems(iInd).controllerGains.decenFSFStabCont{jInd} = K_ijVals{jInd};
                    K{iInd,jInd} = K_ijVals{jInd};
                    K{jInd,iInd} = K_jiVals{jInd};
                end
                
                if ~isStabilizable
                    break
                end
                
            end
            
            if isStabilizable
                Kmat = [];
                for i = 1:1:length(obj.subsystems)
                    Karray = [];
                    for j = 1:1:length(obj.subsystems)
                        Karray = [Karray, K{i,j}];
                    end
                    Kmat = [Kmat; Karray];
                end
                K = Kmat;
            end
            
            % Collect all the coefficients
        end
        
        
        
        
        % Decentralized stable observer design
        function [L,isObserverStable] = decentralizedStableObserverDesign(obj,indexing,solverOptions)
        
            if isempty(indexing)
                indexing = [1:1:length(obj.subsystems)];
            end
            
            for i = 1:1:length(indexing)
                iInd = indexing(i);
                previousSubsystems = indexing(1:i-1);                
                [isObserverStable,L_ii,L_ijVals,L_jiVals] = obj.subsystems(iInd).stableObserverDesign(previousSubsystems, obj.subsystems, solverOptions);
                
                L{iInd,iInd} = L_ii;
                obj.subsystems(iInd).observerGains.decenStabObs{iInd} = L_ii;
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    obj.subsystems(jInd).observerGains.decenStabObs{iInd} = L_jiVals{jInd};
                    obj.subsystems(iInd).observerGains.decenStabObs{jInd} = L_ijVals{jInd};
                    L{iInd,jInd} = L_ijVals{jInd};
                    L{jInd,iInd} = L_jiVals{jInd};
                end
                
                if ~isObserverStable
                    break
                end
                
            end
            
            if isObserverStable
                Lmat = [];
                for i = 1:1:length(obj.subsystems)
                    Larray = [];
                    for j = 1:1:length(obj.subsystems)
                        Larray = [Larray, L{i,j}];
                    end
                    Lmat = [Lmat; Larray];
                end
                L = Lmat;
            end
            
            % Collect all the coefficients
        end
        
        
        % Decentralized DOF stabilization
        function [Ac,Bc,Cc,Dc,isDOFStabilizable] = decentralizedDOFStabilization(obj,indexing,solverOptions)
        
            if isempty(indexing)
                indexing = [1:1:length(obj.subsystems)];
            end
            
            for i = 1:1:length(indexing)
                iInd = indexing(i);
                previousSubsystems = indexing(1:i-1);                
                [isDOFStabilizable,AcVals,BcVals,CcVals,DcVals] = obj.subsystems(iInd).DOFStabilization(previousSubsystems, obj.subsystems, solverOptions);
                
                Ac{iInd,iInd} = AcVals{1}; %A_c,ii
                Bc{iInd,iInd} = BcVals{1};
                Cc{iInd,iInd} = CcVals{1};
                Dc{iInd,iInd} = DcVals{1};
                obj.subsystems(iInd).controllerGains.decenDOFStabContAc{iInd} = AcVals{1};
                obj.subsystems(iInd).controllerGains.decenDOFStabContBc{iInd} = BcVals{1};
                obj.subsystems(iInd).controllerGains.decenDOFStabContCc{iInd} = CcVals{1};
                obj.subsystems(iInd).controllerGains.decenDOFStabContDc{iInd} = DcVals{1};
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    obj.subsystems(iInd).controllerGains.decenDOFStabContAc{jInd} = AcVals{2}{jInd};
                    obj.subsystems(jInd).controllerGains.decenDOFStabContAc{iInd} = AcVals{3}{jInd};                    
                    Ac{iInd,jInd} = AcVals{2}{jInd};
                    Ac{jInd,iInd} = AcVals{3}{jInd};
                    
                    obj.subsystems(iInd).controllerGains.decenDOFStabContBc{jInd} = BcVals{2}{jInd};
                    obj.subsystems(jInd).controllerGains.decenDOFStabContBc{iInd} = BcVals{3}{jInd};                    
                    Bc{iInd,jInd} = BcVals{2}{jInd};
                    Bc{jInd,iInd} = BcVals{3}{jInd};
                    
                    obj.subsystems(iInd).controllerGains.decenDOFStabContCc{jInd} = CcVals{2}{jInd};
                    obj.subsystems(jInd).controllerGains.decenDOFStabContCc{iInd} = CcVals{3}{jInd};                    
                    Cc{iInd,jInd} = CcVals{2}{jInd};
                    Cc{jInd,iInd} = CcVals{3}{jInd};
                    
                    obj.subsystems(iInd).controllerGains.decenDOFStabContDc{jInd} = DcVals{2}{jInd};
                    obj.subsystems(jInd).controllerGains.decenDOFStabContDc{iInd} = DcVals{3}{jInd};                    
                    Dc{iInd,jInd} = DcVals{2}{jInd};
                    Dc{jInd,iInd} = DcVals{3}{jInd};
                end
                
                if ~isDOFStabilizable
                    break
                end
                
            end
            
            if isDOFStabilizable
                Acmat = []; Bcmat = []; Ccmat = []; Dcmat = [];
                for i = 1:1:length(obj.subsystems)
                    AcArray = []; BcArray = []; CcArray = []; DcArray = [];
                    for j = 1:1:length(obj.subsystems)
                        AcArray = [AcArray, Ac{i,j}];
                        BcArray = [BcArray, Bc{i,j}];
                        CcArray = [CcArray, Cc{i,j}];
                        DcArray = [DcArray, Dc{i,j}];
                    end
                    Acmat = [Acmat; AcArray]; 
                    Bcmat = [Bcmat; BcArray];
                    Ccmat = [Ccmat; CcArray];
                    Dcmat = [Dcmat; DcArray];
                end
                Ac = Acmat; Bc = Bcmat; Cc = Ccmat; Dc = Dcmat;
            end
            
            % Collect all the coefficients
        end
        
        
        
        %% Decentralized dissipativity related concepts
        
        
        
        % Check QSR-dissipativity of the networked system
        function output = decentralizedDissipativityAnalysis(obj,dissFrom,dissTo,indexing,solverOptions)
            
            if isempty(indexing)
                indexing = [1:1:length(obj.subsystems)];
            end
            
            for i = 1:1:length(indexing)
                iInd = indexing(i);
                previousSubsystems = indexing(1:i-1);                
                isQSRDissipative = obj.subsystems(iInd).dissipativityAnalysis(dissFrom,dissTo,previousSubsystems, obj.subsystems, solverOptions);
                if ~isQSRDissipative
                    output = false;
                    break
                else
                    output = true;
                end
            end
            
        end
        
        
                
        
        % Dissipating state-feedback controller design
        function [K, isDissipative] = decentralizedFSFDissipativation(obj,indexing,solverOptions)
        
            if isempty(indexing)
                indexing = [1:1:length(obj.subsystems)];
            end
            
            K = [];
            for i = 1:1:length(indexing)
                iInd = indexing(i);
                previousSubsystems = indexing(1:i-1);                
                [isDissipative,K_ii,K_ijVals,K_jiVals] = obj.subsystems(iInd).FSFDissipativation(previousSubsystems, obj.subsystems, solverOptions);
                
                K{iInd,iInd} = K_ii;
                obj.subsystems(iInd).controllerGains.decenFSFDissCont{iInd} = K_ii;
                for j = 1:1:length(previousSubsystems)
                    jInd = previousSubsystems(j);
                    obj.subsystems(jInd).controllerGains.decenFSFDissCont{iInd} = K_jiVals{jInd};
                    obj.subsystems(iInd).controllerGains.decenFSFDissCont{jInd} = K_ijVals{jInd};
                    
                    K{iInd,jInd} = K_ijVals{jInd};
                    K{jInd,iInd} = K_jiVals{jInd};
                end
                
                if ~isDissipative
                    break
                end
            end
            
            if isDissipative
                Kmat = [];
                for i = 1:1:length(obj.subsystems)
                    Karray = [];
                    for j = 1:1:length(obj.subsystems)
                        Karray = [Karray, K{i,j}];
                    end
                    Kmat = [Kmat; Karray];
                end            
                K = Kmat;
            end
                    
        end
        
        
       
    end
end

