%% ----- Homework 5 ----- %%
clc;
clear;
close all;

r1 = 3;
c1 = 3;
r2 = 3;
c2 = 3;
rank_of_A = 3;

% Using the outer product between vectors to control the rank of a matrix.
A = complex(zeros(r1*r2,c1*c2),0);
for i = 1:rank_of_A
    aux1 = randn(r1*r2,1)+ 1i*randn(r1*r2,1); 
    aux2 = randn(1,c1*c2)+ 1i*randn(1,c1*c2); 
    A = A + aux1*aux2; 
end
[U,S,V,rkp] = tensor.KPSVD(A,r1,c1,r2,c2);

% Reconstructing A with full rank.
Ahat = complex(zeros(r1*r2,c1*c2),0);
for r = 1:rkp
    aux1 = U(:,r);
    aux2 = conj(V(:,r));
    Uk{r} = reshape(aux1,[r1 c1]);
    Vk{r}  = reshape(aux2, [r2 c2]); 
    Ahat = Ahat + S(r,r)*tensor.mtx_prod_kron(Uk{r},Vk{r});
end 
disp('Checking the NMSE (dB) between the original matrix A and its'... 
    'reconstruction with KPSVD using full rank:')
Ahat = conj(Ahat);
nmse = (norm(A - Ahat,'fro')^2)/(norm(A,'fro')^2);
nmse = 20*log10(nmse)

% Reconstructing A with deficient rank.
Ahat = complex(zeros(r1*r2,c1*c2),0);
for r = 1:3
    aux1 = U(:,r);
    aux2 = conj(V(:,r));
    Uk{r} = reshape(aux1,[r1 c1]);
    Vk{r}  = reshape(aux2, [r2 c2]); 
    Ahat = Ahat + S(r,r)*tensor.mtx_prod_kron(Uk{r},Vk{r});
end  
disp('Checking the NMSE (dB) between the original matrix A and its'...
    'reconstruction with KPSVD using deficient rank:')
Ahat = conj(Ahat);
nmse = (norm(A - Ahat,'fro')^2)/(norm(A,'fro')^2);
nmse = 20*log10(nmse)
