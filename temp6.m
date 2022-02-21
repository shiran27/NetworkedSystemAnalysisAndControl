% Debugging stabilizing controllers actual vs centralized vs decentralized
clear all
close all 
clc


% Create a network object
network = Network(0);

% Number of subsystems
numOfSubsystems = 4; 
dimentionOfSpace = 2;
sizeOfSpace = 1;
communicationRadius = 0.9;

% subsystem dims
for i = 1:1:numOfSubsystems
    dims{i}.n = 4; % x
    dims{i}.p = 3; % u
    dims{i}.q = 2; % w  
    dims{i}.m = 1; % y
end


errorCount11 = 0;
errorCount12 = 0;
goodCount1 = 0;
errorCount21 = 0;
errorCount22 = 0;
goodCount2 = 0;
numOfTests = 200;
iList1 = [];
iList2 = [];
iList3 = [];
for i = 100 %1:1:numOfTests %48, 782 %544
    i
    % rand('seed',7);
    rng(i)

    % Create the network and plot it
    network = network.loadARandomNetwork(numOfSubsystems,dimentionOfSpace,sizeOfSpace,communicationRadius,dims);
    % network = network.loadTheCustomNetwork();
%     close all
%     network.drawNetwork(1,true);

    [bestIndexing, minCost, worstIndexing, maxCost, basicIndexingCost] = network.findOptimumIndexing()
    % network.drawIndexing(bestIndexing)
    % network.drawIndexing(worstIndexing)

    % Creating the state space representation of the system
    [A,B,C,D,E,F,x] = network.getNetworkMatrices();
    networkedSystem = ss(A,B,C,D);
    isStable1 = isstable(networkedSystem);


    % Check the stability of the networked system in a distributed manner
    isStable2 = network.centralizedStabilityTest();
    [K1, isStabilizable11] = network.designGlobalStabilizingSFBControllers()
    if isStabilizable11
        isStabilizable12 = ~any(real(eig(A+B*K1))>0)
        if isStabilizable11~=isStabilizable12 %actual vs centralized lmi
            if isStabilizable11>isStabilizable12
                errorCount11 = errorCount11 + 1;
            else
                errorCount12 = errorCount12 + 1;
            end
        else
            goodCount1 = goodCount1 + 1;
        end
    end
    
%     isStable3 = network.checkStability([],1) % distributed check stability method 1
    isStable4 = network.checkStability([],2); % distributed check stability method 2
    [K2, isStabilizable21] = network.designLocalStabilizingSFBControllers([]);
    if isStabilizable21
        isStabilizable22 = ~any(real(eig(A+B*K2))>0);
        if isStabilizable21~=isStabilizable22 %actual vs decentralized lmi
            if isStabilizable21>isStabilizable22
                errorCount21 = errorCount21 + 1;
                iList1 = [iList1, i];
%                 break
            else
                errorCount22 = errorCount22 + 1;
            end
        else
            goodCount2 = goodCount2 + 1;
            iList2 = [iList2, i];
        end
    end
    
    if ~isStable1 && ~isStable4 && isStabilizable21 && isStabilizable22 %
%     hard condition tried 5000
%     if ~isStable4 && isStabilizable21 && isStabilizable22 % hard condition
        iList3 = [iList3, i];
%         network.drawNetwork(1,true);
%         i
%         break
    end
        
    
    
    
end
errorCount11
errorCount12
errorCount21
errorCount22