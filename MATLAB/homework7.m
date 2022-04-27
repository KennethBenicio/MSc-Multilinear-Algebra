%% ----- Homework 7 ----- %%
clc;
clear;
close all;

% Full HOSVD
load('Practice_6_hosvd.mat');
[tenS_hat,U_hat] = tensor.HOSVD_full(tenX);
tenX_hat = tensor.n_mod_prod(tenS_hat,U_hat);

% Checking the orthogonality
reshape(tenS_hat(:,:,1),[],1)'*reshape(tenS_hat(:,:,4),[],1)

nmsex = (norm(tensor.unfold(tenX- tenX_hat,1),'fro')^2)/(norm(tensor.unfold(tenX,1),'fro')^2);
nmsex = 20*log10(nmsex)

nmseu1 = (norm(U1 - U_hat{1},'fro')^2)/(norm(U1,'fro')^2);
nmseu2 = (norm(U2 - U_hat{2},'fro')^2)/(norm(U2,'fro')^2);
nmseu3 = (norm(U3 - U_hat{3},'fro')^2)/(norm(U3,'fro')^2);

% Truncated HOSVD and Denoising
X = randn(8,4,10) + 1i*randn(8,4,10);
Y = randn(5,5,5) + 1i*randn(5,5,5);

[S1,U1] = tensor.HOSVD_truncated(X);
multilinear_rank1 = size(S1);
[S2,U2] = tensor.HOSVD_truncated(Y);
multilinear_rank2 = size(S2);

Xhat = tensor.n_mod_prod(S1,U1);
Yhat = tensor.n_mod_prod(S2,U2);

nmsex = (norm(tensor.unfold(X- Xhat,1),'fro')^2)/(norm(tensor.unfold(X,1),'fro')^2);
nmsex = 20*log10(nmsex)
nmsey = (norm(tensor.unfold(Y- Yhat,1),'fro')^2)/(norm(tensor.unfold(Y,1),'fro')^2);
nmsey = 20*log10(nmsey)