
A_ii = network.subsystems(1).A{1}
setlmis([])  % To initialize the LMI description
P = lmivar(1,[size(A_ii,1), 1]) % P is the variable, 1: square symmetric, size(A,1) gives the size and 1 gives that P is a full matrix 

lmiterm([-1, 1, 1, P],-1,A_ii,'s') % defines -PA-A'P; here "s" flag is used to get the symmetric part of PA.
% lmiterm([-1, 1, 2, 0],1)           % defines 0 in M matrix
% lmiterm([-1, 2, 1, 0],1)           % defines 0 in M matrix
lmiterm([-2, 1, 1, P],1,1)
lmisys = getlmis

[tmin,Psol] = feasp(lmisys)
P = dec2mat(lmisys,Psol,P)

% checking
eig(P)
M = -P*A_ii-A_ii'*P
eig(M)