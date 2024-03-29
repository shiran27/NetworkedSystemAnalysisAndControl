classdef Edge < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    properties
        index
        subsystems
        locations
        enabled
    end
    
    methods
        function obj = Edge(index,subsystems,locations)
            obj.index = index;
            obj.subsystems = subsystems;
            obj.locations = locations;
        end
        
        function outputArg = drawEdge(obj)
            hold on
            
            x = obj.locations(:,1);
            y = obj.locations(:,2);
            
            p1 = [x(1) y(1)];
            p2 = [x(2) y(2)];
            dp = p2-p1;
            
            if ~obj.enabled 
%                 plot(obj.locations(:,1),obj.locations(:,2),'k','LineWidth',2, ...
%                     'CreateFcn', @(l, e) set(l, 'Color', [0 0 0 .1]));
            else
%                 plot(obj.locations(:,1),obj.locations(:,2),'k','LineWidth',2);
                q = quiver(p1(1),p1(2),dp(1),dp(2),0,'filled','Color','k','LineWidth',1.25,'MaxHeadSize',0.15/norm(dp));
                
            end
        end
        
        function distanceOutput = getLength(obj)
            distanceOutput = norm(obj.locations(1,:)-obj.locations(2,:),2);
            if ~obj.enabled
                distanceOutput = max(10*distanceOutput,2);
            end
                
        end
        
        
        
        
        
    end
    
end

