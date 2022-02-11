classdef Network < handle
    % This is the class for networks (networked systems)
        
    properties
        index
        subsystems = [] % list of subsystem objects
        edges = []  % list of edge objects
        distanceMatrix
        neighbors = {}
    end
    
    methods
        function obj = Network(index)
            % Construct an instance of this class
            obj.index = index;
        end
        
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
        
        
        function output = loadNeighbors(obj)
            
            for currentSubsystemId = 1:1:length(obj.subsystems) 
                neighborSet = [];
                for subsystemId = 1:1:length(obj.subsystems)
                    % distance will be large if not connected
                    if obj.distanceMatrix(subsystemId,currentSubsystemId) < 1.5  
                        neighborSet = [neighborSet, subsystemId];
                    end
                end
                obj.neighbors{currentSubsystemId} = neighborSet;
                obj.subsystems(currentSubsystemId).neighbors = neighborSet;
            end
        
        end
        
        
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
            height = 0.01

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
        
        function outputArg = finishUpdate(obj)
            for i = 1:1:length(obj.subsystems)
                obj.subsystems(i).finishUpdate();
            end
        end
        
        
        
    end
end

