% Debugging dissipativity tests as well as dissipativity controllers: centralized vs decentralized 
clear all
close all 
clc


% Create a network object
network = Network(0);

% Number of subsystems
numOfSubsystems = 3; 
dimentionOfSpace = 2;
sizeOfSpace = 1;
communicationRadius = 0.9;

% subsystem dims
for i = 1:1:numOfSubsystems
    dims{i}.n = 2; % x
    dims{i}.p = 1; % u
    dims{i}.q = 1; % w  
    dims{i}.m = 1; % y
end


errorCount11 = 0;
errorCount12 = 0;
goodCount1 = 0;
errorCount21 = 0;
errorCount22 = 0;
errorCount31 = 0;
errorCount32 = 0;
errorCount41 = 0;
errorCount42 = 0;
goodCount2 = 0;
numOfTests = 5000;
iList1 = [];
iList2 = [];
iList3 = [];
iList4 = [];
iList5 = [];
for i = 1:1:numOfTests %48, 782 %544
    i
    % rand('seed',7);
    rng(i)

    % Create the network and plot it
    network = network.loadARandomNetwork(numOfSubsystems,dimentionOfSpace,sizeOfSpace,communicationRadius,dims);
    % network = network.loadTheCustomNetwork();
%     close all
%     network.drawNetwork(1,true);

%     [bestIndexing, minCost, worstIndexing, maxCost, basicIndexingCost] = network.findOptimumIndexing()
    % network.drawIndexing(bestIndexing)
    % network.drawIndexing(worstIndexing)

    % Creating the state space representation of the system
    [A,B,C,D,E,F,x] = network.getNetworkMatrices();
    networkedSystem = ss(A,E,C,F);
    isPassive1 = isPassive(networkedSystem);

% Check QSR dissipativity from w to y
% construct Q,S,R here
    [Q,S,R] = network.getSomeQSRMatrices('passive'); % strictly passive, random, passive
    network.storeQSRMatrices(Q,S,R); % store information locally at each sub-system
    isQSRDissipative1 = network.centralizedQSRDissipativityTest();
    isQSRDissipative2 = network.checkQSRDissipativity([]);
    [K1, isDissipative1] = network.designGlobalDissipatingSFBControllers();
    [K2, isDissipative2] = network.designLocalDissipatingSFBControllers([]);

    if isPassive1
        iList1 = [iList1, i];
    end
    
    if isQSRDissipative1~=isQSRDissipative2 %Analysis: centralized lmi vs decentralized lmi
        if isQSRDissipative1>isQSRDissipative2
            errorCount11 = errorCount11 + 1;
            iList2 = [iList2, i];
        else
            errorCount12 = errorCount12 + 1;
        end
    end
    
    if isDissipative1~=isDissipative2 %Synthesis: centralized lmi vs decentralized lmi
        if isDissipative1>isDissipative2
            errorCount21 = errorCount21 + 1;
            iList3 = [iList3, i];
        else
            errorCount22 = errorCount22 + 1;
        end
    end
    
    if isPassive1~=isQSRDissipative1 %Analysis actual passivity vs centralized lmi
        if isPassive1>isQSRDissipative1
            errorCount31 = errorCount31 + 1;
            iList4 = [iList4, i];
        else
            errorCount32 = errorCount32 + 1;
        end
    end
    
    if isPassive1~=isQSRDissipative2 %Analysis: actual vs decentralized lmi
        if isPassive1>isQSRDissipative2
            errorCount41 = errorCount41 + 1;
            iList5 = [iList5, i];
        else
            errorCount42 = errorCount42 + 1;
        end
    end
        
    
    
    
end
% errorCount11
% errorCount12
% errorCount21
% errorCount22