%% ----- Homework 9 ----- %%
clc;
clear;
close all;

ia = 2;
ib = 3;
ic = 4;
ja = 5;
jb = 6;
jc = 7;

rows = [ia ib ic];
columns = [ja jb jc];

A1 = randn(ia,ja);
A2 = randn(ib,jb);
A3 = randn(ic,jc);

aux = tensor.mtx_prod_kron(A1,A2);
X = tensor.mtx_prod_kron(aux,A3);
Ahat = tensor.MLSKronF(X,rows,columns);

aux = tensor.mtx_prod_kron(Ahat{1},Ahat{2});
Xhat = tensor.mtx_prod_kron(aux,Ahat{3});
nmse = (norm(X - Xhat,'fro')^2)/(norm(X,'fro')^2);
nmse = 20*log10(nmse)