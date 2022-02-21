% Generate the system
clear all
close all 
clc

%% Create a network object
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

rng(5) % 5 seems to be the ideal case

% Create the network and plot it
network = network.loadARandomNetwork(numOfSubsystems,dimentionOfSpace,sizeOfSpace,communicationRadius,dims);
% network = network.loadTheCustomNetwork();
network.shiftLocations(-0.05,-0.2);
network.drawNetwork(1,[false,0,1,0.2,0.8]);


% Creating the state space representation of the system and actual stability and passivity
[A,B,C,D,E,F,x] = network.getNetworkMatrices()
networkedSystem = ss(A,E,C,F);
isStable1 = isstable(networkedSystem)   % 0
isPassive1 = isPassive(networkedSystem) % 0


%% Check stability in a distributed manner, and derive stabilizing controller
isStable2 = network.checkStability([],2)                                    % 0
[K1, isStabilizable21] = network.designLocalStabilizingSFBControllers([])   % 1
isStabilizable22 = ~any(real(eig(A+B*K1))>0)                                % 1

%% Check passivity in a distributed manner and derive a passivating controller
[Q,S,R] = network.getSomeQSRMatrices('passive') % strictly passive, random, passive
network.storeQSRMatrices(Q,S,R); % store information locally at each sub-system
isPassive2 = network.checkQSRDissipativity([])  % 0

[K21, isPassive3] = network.designGlobalDissipatingSFBControllers()         % 1 (marginally)
% [K22, isDissipative2] = network.designLocalDissipatingSFBControllers([]);
K21 = K21/1000;
networkedSystem2 = ss(A+B*K21,E,C,F);                                       
isPassive4 =  isPassive(networkedSystem2)                                   % 1


%% 
i = 826
rng(i)
[Q,S,R] = network.getSomeQSRMatrices('random') % strictly passive, random, passive
network.storeQSRMatrices(Q,S,R); % store information locally at each sub-system
isQSRDissipative1 = network.checkQSRDissipativity([])       % 1




%% Find the optimal indexing scheme
[bestIndexing, minCost, worstIndexing, maxCost, basicIndexingCost] = network.findOptimumIndexing()
% network.drawIndexing(bestIndexing)
% network.drawIndexing(worstIndexing)

%% Rounded values
A1 = round(A,3,'significant')
B1 = round(B,3,'significant')
C1 = round(C,3,'significant')
D1 = round(D,3,'significant')
E1 = round(E,3,'significant')
F1 = round(F,3,'significant')

K11 = round(K1,3,'significant')
K211 = round(K21,3,'significant')
Q1 = round(Q,3,'significant')
S1 = round(S,3,'significant')
R1 = round(R,3,'significant')


% Simulate the uncontrolled system (with some qsr passivity level)



% Simulate the stabilized system


% Simulate the passivated system

