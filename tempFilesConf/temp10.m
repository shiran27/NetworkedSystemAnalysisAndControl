close all
clear all
clc

n = 4
p = 3 % inputs
m = 2 % outputs
q = 1 % disturbance input

rsys = rss(n,m,p) % p inputs
rsys2 = rss(n,m,q) % q inputs

A = rsys.A
B = rsys.B
C = rsys.C
D = rsys.D
E = rsys2.B
F = rsys2.D


%% Stabilizing controller design
% P = sdpvar(n,n);
% K = sdpvar(p,n);
% 
% constraints = [P>=0, transpose(A+B*K)*P + P*(A+B*K) <= 0]
% optimize(constraints)
% Pfeasible = value(P)
% Kfeasible = value(K)
% 
% eig(A+B*Kfeasible)


%% Passivating controller design
rho = 2; nu = 2;
Q = -rho*eye(m)
S = 0.5*eye(m,q)
R = -nu*eye(q)

Qhat = C'*Q*C
Shat = C'*S + C'*Q*F
Rhat = F'*Q*F + F'*S + S'*F + R

P = sdpvar(n,n);
K = sdpvar(p,n);

M11 = -transpose(A+B*K)*P - P*(A+B*K) + Qhat
M12 = -P*E+Shat
M21 = transpose(-P*E+Shat)
M22 = Rhat

constraints = [P>=0, [M11,M12;M21,M22]>=0]
r = optimize(constraints)
feasibility = r.problem

P = value(P)
Kfeasible = value(K)

