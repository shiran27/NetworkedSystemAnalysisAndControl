% test the approach

clear all
close all 
clc

rng(6)
numOfSubsystems = 7; 
permList = perms([1:1:numOfSubsystems]);
L = size(permList,1);

numOftests = 100;
errorCount = 0;
positiveTests = 0;
negativeTests = 0;
positiveErrors = 0;
negativeErrors = 0;
halt = false;
for i = 1:1:numOftests

    % Create a network object
    network = Network(0);

    % Number of subsystems
    
    dimentionOfSpace = 2;
    sizeOfSpace = 1;
    communicationRadius = 0.7;

    % Create the network and plot it
    network = network.loadARandomNetwork(numOfSubsystems,dimentionOfSpace,sizeOfSpace,communicationRadius);

%     network.drawNetwork(1);

    
    for j = 1:1:numOftests
        
        bestIndexing = permList(randi(L,1),:); % select a random indexing
        %     network.drawIndexing(bestIndexing)
        
        testNetworkMatrix = network.getARandomNetworkMatrix(4);
        isPositiveDefinite1 = all(eig(testNetworkMatrix)>0);
        isPositiveDefinite2 = network.checkPositiveDefiniteness(testNetworkMatrix,bestIndexing);
        
        if isPositiveDefinite1
            positiveTests = positiveTests + 1;
        else
            negativeTests = negativeTests + 1;
        end
        
        if isPositiveDefinite1~=isPositiveDefinite2
            errorCount = errorCount + 1;        
            disp(['Error at i=',num2str(i),', j=',num2str(j),'.'])
            
            if isPositiveDefinite1 
                positiveErrors = positiveErrors + 1;
            else
                negativeErrors = negativeErrors + 1;
            end
            
%             network.drawNetwork(errorCount);
%             network.drawIndexing(bestIndexing);
%             testNetworkMatrix;
%             halt = true;
%             break

        end
        
    end
    
    if halt
        break
    end
    
end
positiveTests
negativeTests
errorCount
percentageError = errorCount*100/(positiveTests+negativeTests)

positiveErrors
negativeErrors
