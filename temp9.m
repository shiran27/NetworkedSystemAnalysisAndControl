clear all
close all
clc

% Checking custom LMI variables
i = 4
setlmis([])
dim_p=3;
dim_n=2;

% eyeMat = []
% for k = 1:1:(i-1)
%     eyeMatRow = [];
%     for l = 1:1:(i-1)
%         if k==l
%             eyeMatRow = [eyeMatRow,eye(dim_p,dim_n)];
%         else
%             eyeMatRow = [eyeMatRow,zeros(dim_p,dim_n)];
%         end
%     end
%     eyeMat = [eyeMat;eyeMatRow];
% end

% Defining Q_i = diag(Q_i1,I,I,...) with Q_i1 = inv(P_i)
[I_nn,m,sI_nn] = lmivar(1,[dim_n, 0])
[Qi1,m,sQi1] = lmivar(1,[dim_n, 1]);
sQi = [sQi1];
for k = 1:1:(i-1)
    sQi = [sQi, zeros(k*dim_n,dim_n); zeros(dim_n,k*dim_n), sI_nn]; %
end
[Qi,n,sQi] = lmivar(3,sQi);


% Defining L_i = diag(L_ii,L_i1,L_i2,...,L_i,i-1) with L_ii = K_iiP_i^{-1}, L_i1 = K_i1, L_i2 = K_i2,...
[Lij,n,sLij] = lmivar(2,[dim_p,dim_n]); 
sLi = [sLij];
for j = 1:1:(i-1)
    [Lij,~,sLij] = lmivar(2,[dim_p dim_n]);
    Lijvars{j} = Lij
%     eval(['[Li',num2str(j),',~,sLij] = lmivar(2,[dim_p dim_n]);'])
    
    sLi = [sLi, zeros(j*dim_p,dim_n); zeros(dim_p,j*dim_n), sLij];
end
[Li,n,sLi] = lmivar(3,sLi);


% Defining L_ji = diag(L_ji1,I,I,...) with L_ji1 = K_ji*inv(P_i)
[I_pn,m,sI_pn] = lmivar(2,[dim_p,dim_n])
for j = 1:1:(i-1)
    [Lji1,m,sLji1] = lmivar(2,[dim_p dim_n]);
    sLji = [sLji1];
    for k = 1:1:(i-1)
        sLji = [sLji, zeros(k*dim_p,dim_n); zeros(dim_p,k*dim_n), sI_pn];
    end
    % [Lji,n,sLji] = lmivar(3,sLji);
    eval(['[L',num2str(j),'i,n,sL',num2str(j),'i] = lmivar(3,sLji);'])
end

% [Lji,n,sLji] = lmivar(1,[5 1;1 0;3 0]);

M3 = rand(length(sQi));
M3 = M3+M3';
M41 = rand(length(sL1i));
M41 = M41+M41';
Z = zeros(length(sQi),length(sL1i))

% Defining 
lmiterm([-1, 1, 1, Qi],1,1);
lmiterm([-1, 1, 2, -Qi],1,1);
lmiterm([-1, 1, 3, -L1i],1,1);
lmiterm([-1, 2, 1, Qi],1,1);
lmiterm([-1, 3, 1, L1i],1,1);
lmiterm([-1, 2, 2, 0],M3);
lmiterm([-1, 2, 3, 0],Z);
lmiterm([-1, 3, 2, 0],Z');
lmiterm([-1, 3, 3, 0],M41);
lmiterm([-2, 1, 1, Qi],1,1);

lmisys = getlmis;
lmisys = setmvar(lmisys,I_nn,1);
lmisys = setmvar(lmisys,I_pn,eye(dim_p,dim_n));
[tmin,sol] = feasp(lmisys)
Qisol = dec2mat(lmisys,sol,Qi) % for some reason preset variables not appear in Q (they appear as zeros)
Lisol = dec2mat(lmisys,sol,Li)

Li1sol = dec2mat(lmisys,sol,Lijvars{1})
Li2sol = dec2mat(lmisys,sol,Lijvars{2})