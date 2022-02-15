clear all
close all 
clc

% rand('seed',7);
rng(6)

% Create a network object
network = Network(0);

% Number of subsystems
numOfSubsystems = 5; 
dimentionOfSpace = 2;
sizeOfSpace = 1;
communicationRadius = 0.7;

% Create the network and plot it
network = network.loadARandomNetwork(numOfSubsystems,dimentionOfSpace,sizeOfSpace,communicationRadius);
% network = network.loadTheCustomNetwork();

network.drawNetwork(1);

[bestIndexing, minCost, worstIndexing, maxCost] = network.findOptimumIndexing()
network.drawIndexing(bestIndexing)
% network.drawIndexing(worstIndexing)


% Creating the state space representation of the system
[A,B,C,D,E,F,x] = network.getNetworkMatrices();
networkedSystem = ss(A,B,C,D);
isStable1 = isstable(networkedSystem)


% Global controller design to stabilize
% Q = 10000*eye(size(A,1));
% R = 10000*eye(size(B,2));
% [K,~,~] = lqr(A,B,Q,R);
% network.assignLocalControllers(K);
% newA = A-B*K; 
% newA(A==0)=0;
% networkedSystem2 = ss(newA,B,C,D);
% isStable2 = isstable(networkedSystem2)


% Check the stability of the networked system in a distributed manner
isStable = network.centralizedStabilityTest()
isStable3 = network.checkStability(bestIndexing,1)
isStable4 = network.checkStability([],1)
isStable5 = network.checkStability(bestIndexing,2)
isStable6 = network.checkStability([],2)
% network.changeNetworkParameters(newA)
% isStable7 = network.checkStability(bestIndexing,1)
% isStable8 = network.checkStability([],1)
% isStable9 = network.checkStability(bestIndexing,2)
% isStable10 = network.checkStability([],2)


% Check QSR dissipativity from w to y
% construct Q,S,R here
[Q,S,R] = network.getSomeQSRMatrices('strictly passive'); % strictly passive, random, passive
network.storeQSRMatrices(Q,S,R); % store information locally at each sub-system
isQSRDissipative1 = network.centralizedQSRDissipativityTest()
isQSRDissipative2 = network.checkQSRDissipativity(bestIndexing)

% Random network matrix positive definiteness test
% testNetworkMatrix = network.getARandomNetworkMatrix(4);
% isPositiveDefinite1 = all(eig(testNetworkMatrix)>0);
% isPositiveDefinite2 = network.checkPositiveDefiniteness(testNetworkMatrix,bestIndexing);


% Find feedback controls: in a centralized manner
% K = network.designGlobalStabilizingSFBControllers()
% network.assignLocalControllers(-K);
K = network.designGlobalDissipatingSFBControllers()
network.assignLocalControllers(K);

% Find feedback controls: in a decentralized manner
% K = network.designLocalStabilizingSFBControllers(bestIndexing)
% K = network.designLocalDissipatingSFBControllers(bestIndexing)



% Simulating the netwoked system
timeResolution = 0.001;
periodT = 5; % 50, 200, 500, 
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
        time = t;
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