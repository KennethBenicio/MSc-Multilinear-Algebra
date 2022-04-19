%% ----- Homework 5 ----- %%
clc;
clear;
close all;

% My tensor for testing
disp('Generic Tensor:');
tsr_xxx = load('test_unfold_fold.mat').xxx
dimension = size(tsr_xxx);

% Unfolding
disp('1st mode unfolding:');
tsr_xxx_1 = tensor.unfold(tsr_xxx,1)
disp('2nd mode unfolding:');
tsr_xxx_2 = tensor.unfold(tsr_xxx,2)
disp('3rd mode unfolding:');
tsr_xxx_3 = tensor.unfold(tsr_xxx,3)

% Folding
disp('1st mode folding:');
tsr_xxx_1 = tensor.fold(tsr_xxx_1,dimension,1)
disp('2nd mode folding:');
tsr_xxx_2 = tensor.fold(tsr_xxx_2,dimension,2)
disp('3rd mode folding:');
tsr_xxx_3 = tensor.fold(tsr_xxx_3,dimension,3)

% N-mode product
Ytest = tensor.n_mode_prod(X,Matrices,Modes);
nmsey = (norm(Y- Ytest,'fro')^2)/(norm(Y,'fro')^2);
nmsey = 20*log10(nmsey)
