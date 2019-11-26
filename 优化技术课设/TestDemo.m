clear
clc

%Test using two-dimensional constraints
% A=[-1 1;1 1;0 -1;0 1];
% b=[3;3;0;2];
%Test using three-dimensional constraints
A=[0 0 -1;-1 0 1;1 0 1;0 -1 1;0 1 1];
b=[0;3;3;3;3];
[L,U,V]=maxvolrectangle(A,b);