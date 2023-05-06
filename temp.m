clear all
close all
clc

% 

O = zeros(2,2);
I = eye(2);
% 
A = [O, I; O, O];
B = [O; I];
C = [I, O];

% State-space representation 
robot = ss(A,B,C,[])