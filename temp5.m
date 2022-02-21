% Debugging stability check actual vs centralized vs decentralized
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
errorCount21 = 0;
errorCount22 = 0;
numOfTests = 5000;
for i = 1:1:numOfTests
    i
    % rand('seed',7);
    rng(i)

    % Create the network and plot it
    network = network.loadARandomNetwork(numOfSubsystems,dimentionOfSpace,sizeOfSpace,communicationRadius,dims);
    % network = network.loadTheCustomNetwork();
%     close all
%     network.drawNetwork(1,true);

    [bestIndexing, minCost, worstIndexing, maxCost, basicIndexingCost] = network.findOptimumIndexing()
%     network.drawIndexing(bestIndexing)
    % network.drawIndexing(worstIndexing)

    % Creating the state space representation of the system
    [A,B,C,D,E,F,x] = network.getNetworkMatrices();
    networkedSystem = ss(A,B,C,D);
    isStable1 = isstable(networkedSystem);


    % Check the stability of the networked system in a distributed manner
    isStable2 = network.centralizedStabilityTest();
%     isStable3 = network.checkStability([],1) % distributed check stability method 1
    isStable4 = network.checkStability(bestIndexing,2); % distributed check stability method 2
    
    if isStable1~=isStable2 %actual vs centralized lmi
        if isStable1>isStable2
            errorCount11 = errorCount11 + 1;
        else
            errorCount12 = errorCount12 + 1;
        end
    end
    
    if isStable1~=isStable4 %actual vs decentralized lmi
        if isStable1>isStable4
            errorCount21 = errorCount21 + 1;
        else
            errorCount22 = errorCount22 + 1;
        end
    end
end
errorCount11
errorCount12
errorCount21
errorCount22