%% ----- Homework 5 ----- %%
clc;
clear;
close all;

%N = 2;
N = 3;
R = 3;
ia = 4;
ib = 5;
ic = 6;
%dim = [ia ib];
dim = [ia ib ic];

A1 = randn(ia,R) + 1i*randn(ia,R);
A2 = randn(ib,R) + 1i*randn(ib,R);
A3 = randn(ic,R) + 1i*randn(ic,R);
X = tensor.mtx_prod_kr(A1,A2);
X = tensor.mtx_prod_kr(X, A3);

[Ahat] = tensor.MLSKRF(X,N,dim); 
Xhat = tensor.mtx_prod_kr(Ahat{1},Ahat{2});
Xhat = tensor.mtx_prod_kr(Xhat,Ahat{3});
nmse = (norm(X - Xhat,'fro')^2)/(norm(X,'fro')^2);
nmse = 20*log10(nmse)
