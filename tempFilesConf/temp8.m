% generate a system and try different qsr passivity concepts

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
numOfTests = 50000;
iList1 = [];
iList2 = [];
iList3 = [];
iList4 = [];
iList5 = [];

rng(5) % 5 seems to be the ideal case

% Create the network and plot it
network = network.loadARandomNetwork(numOfSubsystems,dimentionOfSpace,sizeOfSpace,communicationRadius,dims);
% network = network.loadTheCustomNetwork();
close all
network.drawNetwork(1,true);

[bestIndexing, minCost, worstIndexing, maxCost, basicIndexingCost] = network.findOptimumIndexing()
network.drawIndexing(bestIndexing)
% network.drawIndexing(worstIndexing)

% Creating the state space representation of the system
[A,B,C,D,E,F,x] = network.getNetworkMatrices();
networkedSystem = ss(A,E,C,F);
isStable1 = isstable(networkedSystem);
isStable2 = network.checkStability([],2);
[K1, isStabilizable21] = network.designLocalStabilizingSFBControllers([]);
isStabilizable22 = ~any(real(eig(A+B*K1))>0);
[Q,S,R] = network.getSomeQSRMatrices('passive'); % strictly passive, random, passive
network.storeQSRMatrices(Q,S,R); % store information locally at each sub-system
isPassive1 = isPassive(networkedSystem);
isPassive2 = network.checkQSRDissipativity([]);
[K21, isDissipative1] = network.designGlobalDissipatingSFBControllers();
[K22, isDissipative2] = network.designLocalDissipatingSFBControllers([]);
networkedSystem2 = ss(A+B*K21,E,C,F);
isPassive3 =  isPassive(networkedSystem2);

for i = 1:1:numOfTests % [826,9652,17549,40414] with random qsr (without positive def. Q and R)
    i
    rng(i)
% Check QSR dissipativity from w to y
% construct Q,S,R here
    [Q,S,R] = network.getSomeQSRMatrices('random'); % strictly passive, random, passive
    network.storeQSRMatrices(Q,S,R); % store information locally at each sub-system
    isQSRDissipative1 = network.centralizedQSRDissipativityTest();
    isQSRDissipative2 = network.checkQSRDissipativity([]);
    [K3, isDissipative1] = network.designGlobalDissipatingSFBControllers();
    [K4, isDissipative2] = network.designLocalDissipatingSFBControllers([]);
    
    if isQSRDissipative2
        i
        disp("Suitable QSR found!!!")
        iList1 = [iList1, i];
    elseif isDissipative2
        i
        disp("Suitable QSR controller found!!!")
        iList2 = [iList2, i];
%         break
    end
end