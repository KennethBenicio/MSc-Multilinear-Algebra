%% ----- Homework 13 ----- %%
clc;
clear;
close all;

% I need to implement the tensor contraction to test the TTSVD.
Ranks = [5 5 5];
X = randn(5,5,5,5);
[G] = tensor.TTSVD(X,Ranks);