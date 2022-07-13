subsystems = network.subsystems;
previousSubsystems = [1,2,3,4];
iInd = 5;
obj = network.subsystems(iInd)

i = length(previousSubsystems)+1;
iInd = obj.index;
disp(['Checking stability at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
A_ii = obj.A{iInd};


% This subsystem has to talk with all the previosSubsystems
% tildeW_ii = [W_ii, W_i; W_i', M_i] > 0 is required where 
% M_i = inv((scriptD_i*scriptA_i)^{-1}*(scriptD_i)*(scriptD_i*scriptA_i)^{-1}')

% LMI Start
setlmis([]);  % To initialize the LMI description
P = lmivar(1,[size(A_ii,1), 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 

% W_ii term
lmiterm([-1, 1, 1, P],-1,A_ii,'s'); % defines -(-PA-A'P)<0; here "s" flag is used to get the symmetric part of PA.

% M_i term
blockSize = obj.dim_n; 
scriptA_i = [];
scriptD_i = [];
for j = 1:1:length(previousSubsystems)
    jInd = previousSubsystems(j);

    % Getting stored info from jInd to create \mathcal{A}_i and \mathcal{D}_i matrices (their j-th columns)
    tildeW_j = subsystems(jInd).dataToBeDistributed.tildeW;

    Z = zeros(blockSize*(i-1-j),blockSize) % (i-1)-j blocks of blockSizeXblockSize zero matrices
    z = zeros(blockSize*(j-1),blockSize);
    if j==1
        tildeW_jj = tildeW_j                    
        scriptA_i = [tildeW_jj; Z];         % The first column of \mathcal{A}_i.
        scriptD_i = [inv(tildeW_jj); Z];    % The first column of \mathcal{D}_i.
    else
        tildeW_jj = tildeW_j(:,blockSize*(j-1)+1:blockSize*j);   % last blockSizeXblockSize block in the row block vector
        tildeW_j  = tildeW_j(:,1:blockSize*(j-1));                % first (j-1) blocks in the row block vector
        scriptA_i = [scriptA_i, [tildeW_j'; tildeW_jj ; Z]];    % The j-th column of \mathcal{A}_i.
        scriptD_i = [scriptD_i, [z; inv(tildeW_jj); Z]];         % The j-th column of \mathcal{D}_i.
    end                    
end
disp(['Data at ',num2str(iInd),' after ',num2str(previousSubsystems),'.'])
scriptD_i
scriptA_i

M1_i = inv(scriptD_i*scriptA_i);
% M_i = inv(M1_i*scriptD_i*M1_i'); % THis fills (i-1)x(i-1) blocks in the LMI
M_i = scriptA_i'*scriptD_i*scriptA_i


% Filling M_i
for x=1:1:length(previousSubsystems)
    m_x = (x-1)*blockSize; % exact index
    n_x = 1+(x-1)+1;           % block index
    for y=1:1:length(previousSubsystems)
        m_y = (y-1)*blockSize;  % exact index
        n_y = 1+(y-1)+1;           % block index
        M_i_xy = M_i(m_x+1:m_x+blockSize, m_y+1:m_y+blockSize);
        lmiterm([-1, n_x, n_y, 0],M_i_xy);
    end
end


% W_i term
for j = 1:1:length(previousSubsystems)
    jInd = previousSubsystems(j);

    % W_ij term
    A_ij = obj.A{jInd};
    A_ji = subsystems(jInd).A{iInd};                    
    P_j = subsystems(jInd).dataToBeDistributed.P;
    n_j = 1+(j-1)+1;

    % Filling W_ij = -A_{ji}^T P_j - P_iA_{ij}
    lmiterm([-1, 1, n_j, P],-1,A_ij);
    matTemp1 = -A_ji'*P_j;
    lmiterm([-1, 1, n_j, 0],matTemp1);

    % Filling W_ij^T = (-A_{ji}^T P_j - P_iA_{ij})^T
    lmiterm([-1, n_j, 1, P],-A_ij',1);
    lmiterm([-1, n_j, 1, 0],matTemp1');
end

% P term
lmiterm([-2, 1, 1, P],1,1); % defines -P<0

% Completing the LMI
lmisys = getlmis;
[tmin,Psol] = feasp(lmisys); % Solve the LMI system
isFeasible = tmin <= 0; % strictly feasible if this is satisfied
P_i = dec2mat(lmisys,Psol,P); % This needs to be stored
% LMI End

% Need to compute \tilede{W}_i and \tilde{W}_{ii} for storage

% W_i term
W_i = [];
for j = 1:1:length(previousSubsystems)
    jInd = previousSubsystems(j);

    % W_ij term
    A_ij = obj.A{jInd};
    A_ji = subsystems(jInd).A{iInd};
    P_j = subsystems(jInd).dataToBeDistributed.P;

    % Filling W_ij = -A_{ji}^T P_j - P_iA_{ij}
    W_ij = -P_i*A_ij - A_ji'*P_j;
    W_i = [W_i, W_ij];
end
tildeW_i = W_i*M1_i;

% W_ii term: W_ii = -A_{ii}^T P_i - P_iA_{ii}
W_ii = -A_ii'*P_i - P_i*A_ii;

tildeW_ii = W_ii - tildeW_i*scriptD_i*tildeW_i'; % Note that here, \tilde{W}_ii, W_ii, \tilde{W}_i are different.
tildeW_i = [tildeW_i, tildeW_ii];

obj.dataToBeDistributed.tildeW = tildeW_i; % Storing
obj.dataToBeDistributed.P = P_i; % Storing
if ~isFeasible
    disp(['LMI is not feasible at: ',num2str(iInd),' after ',num2str(previousSubsystems),'.']);
    output = false;
else
    output = true;
end