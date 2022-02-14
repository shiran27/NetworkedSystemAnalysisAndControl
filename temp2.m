obj = network.subsystems(1)
iInd = 1;

A_ii = obj.A{iInd};
C_i = obj.C{iInd};
E_ii = obj.E{iInd};
F_i = obj.F{iInd};
Q_ii = obj.dataToBeDistributed.Q{iInd};
S_ii = obj.dataToBeDistributed.Q{iInd};
R_ii = obj.dataToBeDistributed.Q{iInd};


% LMI Start
setlmis([]);  % To initialize the LMI description
P = lmivar(1,[size(A_ii,1), 1]); % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 

% W_ii1 = -A_{ii}^T P_i - P_iA_{ii} + C_i^\T Q_{ii} C_i
lmiterm([-1, 1, 1, P],-1,A_ii,'s'); 
matTemp = C_i'*Q_ii*C_i;
lmiterm([-1, 1, 1, 0],matTemp);
% W_ii2 = -P_i E_{ii} + C_i^\T S_{ii} + C_i^\T Q_{ii} F_i;
lmiterm([-1, 1, 2, P],-1,E_ii); 
matTemp = C_i'*S_ii + C_i'*Q_ii*F_i;
lmiterm([-1, 1, 2, 0],matTemp); 
% W_ii3 = -E_{ii}^\T P_i + S_{ii}^\T C_i + F_i^\T Q_{ii}^\T C_i;
lmiterm([-1, 2, 1, P],-E_ii',1);
lmiterm([-1, 2, 1, 0],matTemp'); 
% W_ii4 = F_i^\T Q_{ii} F_i + (F_i^\T S_{ii} + S_{ii}^\T F_i) + R_{ii};
matTemp = F_i'*Q_ii*F_i + (F_i'*S_ii+S_ii'*F_i) + R_ii;
lmiterm([-1, 2, 2, 0],matTemp);

% P>0
lmiterm([-2, 1, 1, P],1,1); % defines -P<0
lmisys = getlmis;

[tmin,Psol] = feasp(lmisys); % Solve the LMI system
isFeasible = tmin <= 0; % strictly feasible if this is satisfied
P_i = dec2mat(lmisys,Psol,P); % This needs to be stored