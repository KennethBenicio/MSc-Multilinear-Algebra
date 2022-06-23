%% ----- Homework 4 ----- %%
clc;
clear;
close all;

A = randn(4,2) + 1j*randn(4,2);
B = randn(6,3) + 1j*randn(6,3);
X = tensor.mtx_prod_kron(A,B);
[Ahat,Bhat] = tensor.LSKronF(X,4,2,6,3);
Xhat = tensor.mtx_prod_kron(Ahat,Bhat);

disp('Checking the NMSE (dB) between the original matrix X and its reconstruction with LSKronF:')
nmsex = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
nmsex = 20*log10(nmsex)
disp('Checking the NMSE (dB) between the original matrix A and its estimation:')
nmsea = (norm(A- Ahat,'fro')^2)/(norm(A,'fro')^2);
nmsea = 20*log10(nmsea)
disp('Checking the NMSE (dB) between the original matrix B and its estimation:')
nmseb = (norm(B- Bhat,'fro')^2)/(norm(B,'fro')^2);
nmseb = 20*log10(nmseb)

I = 2;
J = 4;
P = 3;
Q = 5;
SNR = [0 5 10 15 20 25 30];
nmse = zeros(length(SNR),1);
for snr = 1:length(SNR)
    for mc = 1:1000
        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise/2)*(randn(I*J,P*Q) + 1j*randn(I*J,P*Q));
        
        A = randn(I,P) + 1j*randn(I,P);
        B = randn(J,Q) + 1j*randn(J,Q);
        X = tensor.mtx_prod_kron(A,B);
        X_noisy = X + noise;
        
        [Ahat,Bhat] = tensor.LSKronF(X_noisy,I,P,J,Q);
        Xhat = tensor.mtx_prod_kron(Ahat,Bhat);
        aux = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse(snr,1) = nmse(snr,1) + 20*log10(aux);
    end
end
nmse  = nmse/1000;

figure
txt = ['I = ' num2str(I), ', J = ' num2str(J), ', P = ' num2str(P), ' and Q = ' num2str(Q)];
plot(SNR,nmse,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;

I = 4;
J = 8;
P = 3;
Q = 5;
SNR = [0 5 10 15 20 25 30];
nmse = zeros(length(SNR),1);
for snr = 1:length(SNR)
    for mc = 1:1000
        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise/2)*(randn(I*J,P*Q) + 1j*randn(I*J,P*Q));
        
        A = randn(I,P) + 1j*randn(I,P);
        B = randn(J,Q) + 1j*randn(J,Q);
        X = tensor.mtx_prod_kron(A,B);
        X = X + noise;
        
        [Ahat,Bhat] = tensor.LSKronF(X,I,P,J,Q);
        Xhat = tensor.mtx_prod_kron(Ahat,Bhat);
        aux = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse(snr,1) = nmse(snr,1) + 20*log10(aux);
    end
end
nmse  = nmse/1000;

txt = ['I = ' num2str(I), ', J = ' num2str(J), ', P = ' num2str(P), ' and Q = ' num2str(Q)];
plot(SNR,nmse,'-d','color', [0.4660 0.6740 0.1880], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title(['LSKronF performance under imperfect scenario'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
grid on;
saveas(gcf,'hw4.png')