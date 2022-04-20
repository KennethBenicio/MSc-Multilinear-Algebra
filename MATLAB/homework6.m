%% ----- Homework 6 ----- %%
clc;
clear;
close all;

% Full HOSVD
[X] = load('hosvd_test.mat');

[S,U] = HOSVD_full(X);
Xhat = tensor.n_mod_prod(S,U);
nmsex = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
nmsex = 20*log10(nmsex)

% Truncated HOSVD and Denoising
[X,Y] = load('hosvd_denoising.mat');

[S1,U1] = tensor.HOSVD_truncated(X);
multilinear_rank1 = size(S1);
[S2,U2] = tensor.HOSVD_truncated(Y);
multilinear_rank2 = size(S2);

Xhat = tensor.n_mod_prod(S1,U1);
Yhat =;tensor.n_mod_prod(S2,U2);

nmsex = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
nmsex = 20*log10(nmsex)
nmsey = (norm(Y- Yhat,'fro')^2)/(norm(Y,'fro')^2);
nmsey = 20*log10(nmsey)