%% ----- Homework 9 ----- %%
clc;
clear;
close all;

load('homework9_MLSKRF.mat')

N = 3;
dim = [5 4 8];
[Matrices] = tensor.MLSKRF(X,N,dim); 
Xhat = tensor.mtx_prod_kr(Matrices{1},Matrices{2});
Xhat = tensor.mtx_prod_kr(Xhat,Matrices{3});

disp('Checking the NMSE (dB) between the original matrix X and its reconstruction with MLSKRF:')
nmsex = (norm(X - Xhat,'fro')^2)/(norm(X,'fro')^2);
nmsex = 20*log10(nmsex)
disp('Checking the NMSE (dB) between the original matrix A and its estimation:')
nmsea = (norm(A - Matrices{1},'fro')^2)/(norm(A,'fro')^2);
nmsea = 20*log10(nmsea)
disp('Checking the NMSE (dB) between the original matrix B and its estimation:')
nmseb = (norm(B - Matrices{2},'fro')^2)/(norm(B,'fro')^2);
nmseb = 20*log10(nmseb)
disp('Checking the NMSE (dB) between the original matrix C and its estimation:')
nmsec = (norm(C - Matrices{3},'fro')^2)/(norm(C,'fro')^2);
nmsec = 20*log10(nmsec)

N = 3;
R  = 4;
I = 2;
J = 3;
K = 4;
dim = [I J K];
SNR = [0 5 10 15 20 25 30];
nmse = zeros(length(SNR),1);
for snr = 1:length(SNR)
    for mc = 1:1000
        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise/2)*(randn(I*J*K,R) + 1j*randn(I*J*K,R));
        
        A = randn(I,R) + 1j*randn(I,R);
        B = randn(J,R) + 1j*randn(J,R);
        C = randn(K,R) + 1j*randn(K,R);
        X = tensor.mtx_prod_kr(A,B);
        X = tensor.mtx_prod_kr(X,C);
        X_noisy = X + noise;
        
        [Matrices] = tensor.MLSKRF(X_noisy,N,dim);
        Xhat = tensor.mtx_prod_kr(Matrices{1},Matrices{2});
        Xhat = tensor.mtx_prod_kr(Xhat,Matrices{3});
        aux = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse(snr,1) = nmse(snr,1) + 20*log10(aux);
    end
end
nmse  = nmse/1000;

figure
txt = ['I1 = ' num2str(I), ', I2 = ' num2str(J), ', I3 = ' num2str(K), ' and R = ' num2str(R)];
plot(SNR,nmse,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
title(['MLSKRF performance under imperfect scenario'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
grid on;
saveas(gcf,'hw9.png')
