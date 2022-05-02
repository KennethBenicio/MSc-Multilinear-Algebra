%% ----- Homework 3 ----- %%
clc;
clear;
close all;

A = randn(4,2) + 1j*randn(4,2);
B = randn(6,2) + 1j*randn(6,2);
X = tensor.mtx_prod_kr(A,B);
[Ahat,Bhat] = tensor.LSKRF(X,4,6);
Xhat = tensor.mtx_prod_kr(Ahat,Bhat);

disp('Checking the NMSE (dB) between the original matrix X and its reconstruction with LSKRF:')
nmsex = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
nmsex = 20*log10(nmsex)
disp('Checking the NMSE (dB) between the original matrix A and its estimation:')
nmsea = (norm(A- Ahat,'fro')^2)/(norm(A,'fro')^2);
nmsea = 20*log10(nmsea)
disp('Checking the NMSE (dB) between the original matrix B and its estimation:')
nmseb = (norm(B- Bhat,'fro')^2)/(norm(B,'fro')^2);
nmseb = 20*log10(nmseb)

I = 10;
J = 10;
R =  4;
SNR = [0 5 10 15 20 25 30];
nmse = zeros(length(SNR),1);
for snr = 1:length(SNR)
    for mc = 1:1000
        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise/2)*(randn(I*J,R) + 1j*randn(I*J,R));
        
        A = randn(I,R) + 1j*randn(I,R);
        B = randn(J,R) + 1j*randn(J,R);
        X = tensor.mtx_prod_kr(A,B);
        X_noisy = X + noise;
        
        [Ahat,Bhat] = tensor.LSKRF(X_noisy,I,J);
        Xhat = tensor.mtx_prod_kr(Ahat,Bhat);
        aux = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse(snr,1) = nmse(snr,1) + 20*log10(aux);
    end
end
nmse  = nmse/1000;

figure
txt = ['I = ' num2str(I), ', J = ' num2str(J), ' and R = ' num2str(R)];
plot(SNR,nmse,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
title(['LSKRF performance under imperfect scenario'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
grid on;
saveas(gcf,'hw3a1.png')

I = 30;
J = 10;
R =  4;
SNR = [0 5 10 15 20 25 30];
nmse = zeros(length(SNR),1);
for snr = 1:length(SNR)
    for mc = 1:1000
        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise/2)*(randn(I*J,R) + 1j*randn(I*J,R));
        
        A = randn(I,R) + 1j*randn(I,R);
        B = randn(J,R) + 1j*randn(J,R);
        X = tensor.mtx_prod_kr(A,B);
        X = X + noise;
        
        [Ahat,Bhat] = tensor.LSKRF(X,I,J);
        Xhat = tensor.mtx_prod_kr(Ahat,Bhat);
        aux = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse(snr,1) = nmse(snr,1) + 20*log10(aux);
    end
end
nmse  = nmse/1000;

figure
txt = ['I = ' num2str(I), ', J = ' num2str(J), ' and R = ' num2str(R)];
plot(SNR,nmse,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
title(['LSKRF performance under imperfect scenario'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
grid on;
saveas(gcf,'hw3a2.png')