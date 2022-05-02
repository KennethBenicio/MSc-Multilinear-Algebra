%% ----- Homework 9 ----- %%
clc;
clear;
close all;

load('homework9_MLSKronF.mat')

rows = [4 4 6];
columns = [3 2 5];
[Matrices] = tensor.MLSKronF(X,rows,columns);

aux = tensor.mtx_prod_kron(Matrices{1},Matrices{2});
Xhat = tensor.mtx_prod_kron(aux,Matrices{3});
disp('Checking the NMSE (dB) between the original matrix X and its reconstruction with MLSKRF:')
nmse = (norm(X - Xhat,'fro')^2)/(norm(X,'fro')^2);
nmse = 20*log10(nmse)
disp('Checking the NMSE (dB) between the original matrix A and its estimation:')
nmsea = (norm(A - Matrices{1},'fro')^2)/(norm(A,'fro')^2);
nmsea = 20*log10(nmsea)
disp('Checking the NMSE (dB) between the original matrix B and its estimation:')
nmseb = (norm(B - Matrices{2},'fro')^2)/(norm(B,'fro')^2);
nmseb = 20*log10(nmseb)
disp('Checking the NMSE (dB) between the original matrix C and its estimation:')
nmsec = (norm(C - Matrices{3},'fro')^2)/(norm(C,'fro')^2);
nmsec = 20*log10(nmsec)

ia = 2;
ib = 3;
ic = 4;
ja = 2;
jb = 3;
jc = 4;
rows = [ia ib ic];
columns = [ja jb jc];
SNR = [0 5 10 15 20 25 30];
nmse = zeros(length(SNR),1);
for snr = 1:length(SNR)
    for mc = 1:1000
        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise/2)*(randn(ia*ib*ic,ja*jb*jc) + 1j*randn(ia*ib*ic,ja*jb*jc));
        
        A = randn(ia,ja) + 1i*randn(ia,ja);
        B = randn(ib,jb) + 1i*randn(ib,jb);
        C = randn(ic,jc) + 1i*randn(ic,jc);
        X = tensor.mtx_prod_kron(A,B);
        X = tensor.mtx_prod_kron(X,C);
        X_noisy = X + noise;
        
        [Matrices] = tensor.MLSKronF(X_noisy,rows,columns);
        aux = tensor.mtx_prod_kron(Matrices{1},Matrices{2});
        Xhat = tensor.mtx_prod_kron(aux,Matrices{3});
        aux = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse(snr,1) = nmse(snr,1) + 20*log10(aux);
    end
end
nmse  = nmse/1000;

figure
txt = ['MLSKronF with N = 3'];
plot(SNR,nmse,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
title(['MLSKronF performance under imperfect scenario'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
grid on;
saveas(gcf,'hw9a1.png')