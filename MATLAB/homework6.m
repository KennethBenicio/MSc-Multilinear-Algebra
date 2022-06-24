%% ----- Homework 6 ----- %%
clc;
clear;
close all;

% Tensor for testing
disp('Generic Tensor:');
load('homework6_unfolding_folding.mat')
dimension = size(tenX);

% Kolda example for testing
% load('homework6_kolda_example.mat')
% tenX = xxx;
% dimension = size(tenX);

% Unfolding
disp('X1./unfold(tenX,1):');
tenX_1 = tensor.unfold(tenX,1);
X1./tensor.unfold(tenX,1)
disp('X2./unfold(tenX,2):');
tenX_2 = tensor.unfold(tenX,2);
X2./tensor.unfold(tenX,2)
disp('X3./unfold(tenX,3):');
tenX_3 = tensor.unfold(tenX,3);
X3./tensor.unfold(tenX,3)

% Folding
disp('tenX./fold(unfold(tenX,1),1):');
tenX./tensor.fold(tenX_1,dimension,1)
disp('tenX./fold(unfold(tenX,2),2):');
tenX./tensor.fold(tenX_2,dimension,2)
disp('tenX./fold(unfold(tenX,3),3):');
tenX./tensor.fold(tenX_3,dimension,3)

% N-mode product
load('homework6_n_mode.mat')
tenY_test = tensor.n_mod_prod(tenX,{Z},[1]);
disp('Checking the NMSE (dB) between the original tensor Y and the one'... 
    'after the N-mode product:')
nmsey = (norm(tensor.unfold(tenY- tenY_test,1),'fro')^2)...
    /(norm(tensor.unfold(tenY,1),'fro')^2);
nmsey = 20*log10(nmsey)
