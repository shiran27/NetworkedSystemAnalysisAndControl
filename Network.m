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
        
        leastNeighbors = []     % \Gamma_i 's
        leastOutNeighbors = []  % \Delta_i 's
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
                obj.subsystems(i).loadParameters(obj.subsystems);
%                 obj.subsystems(i).loadStableParameters(obj.subsystems);
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
       
    end
end

