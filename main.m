clear all
close all 
clc

rand('seed',7);

% Create a network object
network = Network(0);

% Number of subsystems
numOfSubsystems = 7;
dimentionOfSpace = 2;
sizeOfSpace = 1;
communicationRadius = 0.7;

% Create the network and plot it
network = network.loadARandomNetwork(numOfSubsystems,dimentionOfSpace,sizeOfSpace,communicationRadius);
network.drawNetwork(1);

[bestIndexing, minCost, worstIndexing, maxCost] = network.findOptimumIndexing()
network.drawIndexing(bestIndexing)
% network.drawIndexing(worstIndexing)


% [A,B,C,D,E,F,x] = network.getNetworkMatrices()




timeResolution = 0.001;
periodT = 1; % 50, 200, 500, 
textHandle1 = text();  textHandle2 = text(); textHandle3 = text(); 
plotMode = true; 
videoMode = false;
data = zeros(periodT/timeResolution+1, numOfSubsystems, 4);
timeSteps = 1;
frameCount = 1;
xCost = 0;
uCost = 0;
for t = 0:timeResolution:periodT
    
    % target update
    for i = 1:1:numOfSubsystems
        [xCostVal, uCostVal, dataVal] = network.subsystems(i).update(timeResolution,network.subsystems); 
        xCost = xCost + xCostVal;
        uCost = uCost + uCostVal;
        
        if plotMode
            data(timeSteps,i,:) = dataVal; %[||x|| ; ||u|| ; ||w|| ; ||y||]
        end
    end
    network.finishUpdate(); % finish updating the state variables in subsystems
    
    
    
    % display & video
    if rem(t,0.1)==0 
        time = t
        if plotMode 
            network.updateNetwork(1,t);
            if t~=0
                delete(textHandle1);
                textHandle1 = text(0.05,0.96,['State Cost: ',num2str(xCost/t,5)],'Color','k','FontSize',10);
            end
            delete(textHandle2);
            textHandle2 = text(0.05,0.92,['Ctrl. Cost: ',num2str(uCost/t,5)],'Color','k','FontSize',10);
            
            delete(textHandle3);
            textHandle3 = text(0.05,0.88,['Time: ',num2str(t,5)],'Color','k','FontSize',10);
            
        end
        
        if videoMode
            frameArray(frameCount) = getframe(gcf);
            frameCount = frameCount + 1;
        end
        
        
    end
    
    timeSteps = timeSteps + 1;
end

if plotMode
    run('graphicsOfSimulation.m')
end