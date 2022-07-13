
clear all
close all
clc

data = [];
for seedVal = 100:1:1000
    
    rng(seedVal) 

    %% Create a network object
    network = Network(0);

    % Number of subsystems
    numOfSubsystems = 5; % 5
    dimentionOfSpace = 2;
    sizeOfSpace = 1;
    communicationRadius = 0.8; %0.9

    % subsystem dims
    for i = 1:1:numOfSubsystems
        dims{i}.n = 2; % x
        dims{i}.p = 1; % u
        dims{i}.q = 1; % w  
        dims{i}.m = 1; % y
        dims{i}.l = 1; % z 
    end

    % Diagonal matrices out of Network Matrices: A,B,C,D,E,F,G,H,J,Q,S,R 
    % (1 if true, -1 % if zero and 0 if false)
    diags = [0,1,1,-1,0,1,1,1,1]; % respectively for 

    % Create the network and plot it
    network = network.loadARandomNetwork(numOfSubsystems,dimentionOfSpace,sizeOfSpace,communicationRadius,dims,diags);
    network.shiftLocations();
    close all
%     network.drawNetwork(1,false);

    % Get network matrices
    [A,B,C,D,E,F,G,H,J,x] = network.getNetworkMatrices();
    n = size(A,1);
    p = size(B,2);
    q = size(E,2);


    % Generate Q,S,R matrices
    dissFrom = 'w';    % Options: 'u', 'w'
    dissTo = 'y';      % Options: 'y', 'z'
    dissType = 'strictly passive';     % Options: strictly passive([IFP, OFP]), passive, L2G(gamma), random 
    dissArgs = [0.2,0.2]; %[0,0.01]
    % dissType = 'passive';
    % dissArgs = 0;
    [Q,S,R] = network.getSomeQSRMatrices(dissFrom,dissTo,dissType,dissArgs); 
    network.storeQSRMatrices(Q,S,R); % Store information locally at each sub-system
    % % %     Central Stability Based Concepts
    % addpath('C:\Program Files\Mosek\9.3\toolbox\R2015a')
    solverOptions = sdpsettings('solver','mosek','verbose',0);
    rng(seedVal);

    % % % % Stability Analysis
    is1 = network.centralizedStabilityAnalysis(solverOptions);

    % % % % Stabilization using FSF control
    [KStab,is2] = network.centralizedFSFStabilization(solverOptions);

    % % % % Stable observer design
    [LStab,is3] = network.centralizedStableObserverDesign(solverOptions);

    % % % % Stabilization using DOF control
    [AcStab,BcStab,CcStab,DcStab,is4] = network.centralizedDOFStabilization(solverOptions);
    
    % % % %     Centralized Dissipativity Based Results
    dissFrom = 'w'; dissTo = 'y';  
    % dissType = 'strictly passive'; %dissArgs = [0.00001,1];
    [Q,S,R] = network.getSomeQSRMatrices(dissFrom,dissTo,dissType,dissArgs); % dissFrom,dissTo,dissType,dissArgs

    % % % % Dissipativity Analysis
    is5 = network.centralizedDissipativityAnalysis(dissFrom,dissTo,Q,S,R,solverOptions);

    % % % % Dissipativation using FSF control
    [KDiss,is6] = network.centralizedFSFDissipativation(dissFrom,dissTo,Q,S,R,solverOptions);

    % % % % Dissipative observer design
    dissFrom = 'w'; dissTo = 'z';  
    % dissType = 'strictly passive'; %dissArgs = [0.00001,1];
    [Q,S,R] = network.getSomeQSRMatrices(dissFrom,dissTo,dissType,dissArgs); % dissFrom,dissTo,dissType,dissArgs
    [LDiss,is7] = network.centralizedDissipativeObserverDesign(dissFrom,dissTo,Q,S,R,solverOptions);

    % % % % Stabilization using DOF control
    dissFrom = 'w'; dissTo = 'z';  
    % dissType = 'strictly passive'; %dissArgs = [0.00001,1];
    [Q,S,R] = network.getSomeQSRMatrices(dissFrom,dissTo,dissType,dissArgs); % dissFrom,dissTo,dissType,dissArgs
    [AcDiss,BcDiss,CcDiss,DcDiss,is8] = network.centralizedDOFDissipativation(dissFrom,dissTo,Q,S,R,solverOptions);
    
    
    
    
    
    % % Decentralized Stability Based Concepts
    % addpath('C:\Program Files\Mosek\9.3\toolbox\R2015a')
    solverOptions = sdpsettings('solver','mosek','verbose',0);
    indexingScheme = []; % use the default scheme
    rng(seedVal);

    % % % % Stability Analysis
    is9 = network.decentralizedStabilityAnalysis([],solverOptions);

    % % % % Stabilization using FSF control
    [KStab,is10] = network.decentralizedFSFStabilization([],solverOptions);

    % % % % Stable observer design
    [LStab,is11] = network.decentralizedStableObserverDesign([],solverOptions);

    % % % % Stabilization using DOF control
    [AcStab,BcStab,CcStab,DcStab,is12] = network.decentralizedDOFStabilization([],solverOptions); 
    % % Decentralized Dissipativity Based Concepts
    dissFrom = 'w'; dissTo = 'y';  
    % dissType = 'strictly passive'; dissArgs = [0.2,0.2];
    [Q,S,R] = network.getSomeQSRMatrices(dissFrom,dissTo,dissType,dissArgs); % dissFrom,dissTo,dissType,dissArgs
    network.storeQSRMatrices(Q,S,R);

    % % % % Dissipativity Analysis
    is13 = network.decentralizedDissipativityAnalysis(dissFrom,dissTo,[],solverOptions);

    % % % % Dissipativation using FSF control
    [KDiss,is14] = network.decentralizedFSFDissipativation(dissFrom,dissTo,[],solverOptions);

    % % % % Dissipative observer design
    dissFrom = 'w'; dissTo = 'z';  
    % dissType = 'strictly passive'; %dissArgs = [0.00001,1];
    [Q,S,R] = network.getSomeQSRMatrices(dissFrom,dissTo,dissType,dissArgs); % dissFrom,dissTo,dissType,dissArgs
    network.storeQSRMatrices(Q,S,R);
    [LDiss,is15] = network.decentralizedDissipativeObserverDesign(dissFrom,dissTo,[],solverOptions);

    % % % % Stabilization using DOF control
    dissFrom = 'w'; dissTo = 'z';  
    % dissType = 'strictly passive'; %dissArgs = [0.00001,1];
    [Q,S,R] = network.getSomeQSRMatrices(dissFrom,dissTo,dissType,dissArgs) % dissFrom,dissTo,dissType,dissArgs
    network.storeQSRMatrices(Q,S,R);
    [AcDiss,BcDiss,CcDiss,DcDiss,is16] = network.decentralizedDOFDissipativation(dissFrom,dissTo,[],solverOptions); 
    
    
    
    dataTemp = [seedVal,is1,is2,is3,is4,is5,is6,is7,is8,seedVal,is9,is10,is11,is12,is13,is14,is15,is16]
    data = [data;dataTemp];
end

data
