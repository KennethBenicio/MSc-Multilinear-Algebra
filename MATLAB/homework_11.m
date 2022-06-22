%% ----- Homework 11 ----- %%
clc;
clear;
close all;

% Random tensor example
% I = 8;
% J = 4;
% K = 5;
% R = 3;
% A = randn(I,R) + 1j*randn(I,R);
% B = randn(J,R) + 1j*randn(J,R);
% C = randn(K,R) + 1j*randn(K,R);
% X = tensor.fold(A*(tensor.mtx_prod_kr(C,B).'),[I J K],1);

load('homework_11_CPD.mat')
X = tenX;
R = 3;

[ia,ib,ic] = size(X);
mode_1 = tensor.unfold(X,1);
mode_2 = tensor.unfold(X,2);
mode_3 = tensor.unfold(X,3);

Ahat = randn(ia,R) + 1j*randn(ia,R);
Bhat = randn(ib,R) + 1j*randn(ib,R);
Chat = randn(ic,R) + 1j*randn(ic,R);

aux = 10000;
error = zeros(1,aux);
error(1) = ((norm((mode_1 - Ahat*(tensor.mtx_prod_kr(Chat,Bhat).')),'fro'))^2)/((norm(mode_1,'fro')^2));
for i = 2:aux
    Bhat = mode_2*pinv((tensor.mtx_prod_kr(Chat,Ahat)).');
    Chat = mode_3*pinv((tensor.mtx_prod_kr(Bhat,Ahat)).');
    Ahat = mode_1*pinv((tensor.mtx_prod_kr(Chat,Bhat)).');
    error(i) = ((norm((mode_1 - Ahat*(tensor.mtx_prod_kr(Chat,Bhat).')),'fro'))^2)/((norm(mode_1,'fro')^2));
    if abs(error(i) - error(i-1)) < eps
        error = error(1:i);
        break;
    else
        continue;
    end
end

disp('Checking the NMSE (dB) between the original tensor X and its reconstruction with MLSKRF:')
Xhat = tensor.fold(Ahat*(tensor.mtx_prod_kr(Chat,Bhat).'),[ia ib ic],1);
nmsex = (norm(tensor.unfold(X- Xhat,1),'fro')^2)/(norm(tensor.unfold(X,1),'fro')^2);
nmsex = 20*log10(nmsex)
disp('Checking the NMSE (dB) between the original matrix A and its estimation:')
nmsea = (norm(A - Ahat,'fro')^2)/(norm(A,'fro')^2);
nmsea = 20*log10(nmsea)
disp('Checking the NMSE (dB) between the original matrix B and its estimation:')
nmseb = (norm(B - Bhat,'fro')^2)/(norm(B,'fro')^2);
nmseb = 20*log10(nmseb)
disp('Checking the NMSE (dB) between the original matrix C and its estimation:')
nmsec = (norm(C - Chat,'fro')^2)/(norm(C,'fro')^2);
nmsec = 20*log10(nmsec)

figure
txt = ['\bf Convergence between iterations'];
plot(1:i,20*log10(error),'-','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
title(['ALS convergence: I = ' num2str(ia), ', J = ' num2str(ib), ', K = ' num2str(ic), ' and R = ' num2str(R)])
xlabel('Iterations')
ylabel('Error (dB)')
legend_copy = legend("location", "southwest");
set(legend_copy,'Interpreter','tex','location','southwest',"fontsize", 12)
grid on;
saveas(gcf,'hw11a1.png')

%% ALS in the presence of noisy signals

I = 8;
J = 4;
K = 5;
R = 3;
SNR = [0 5 10 15 20 25 30];
nmse1 = zeros(length(SNR),1);
nmse2 = zeros(length(SNR),1);
nmse3 = zeros(length(SNR),1);
nmse4 = zeros(length(SNR),1);
for snr = 1:length(SNR)
    for mc = 1:1000
        A = randn(I,R) + 1j*randn(I,R);
        B = randn(J,R) + 1j*randn(J,R);
        C = randn(K,R) + 1j*randn(K,R);
        X = tensor.fold(A*(tensor.mtx_prod_kr(C,B).'),[I J K],1);

        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise/2)*(randn(I,J,K) + 1j*randn(I,J,K));
        
        X = X + noise;
        
        [ia,ib,ic] = size(X);
        mode_1 = tensor.unfold(X,1);
        mode_2 = tensor.unfold(X,2);
        mode_3 = tensor.unfold(X,3);
        Ahat = randn(ia,R) + 1j*randn(ia,R);
        Bhat = randn(ib,R) + 1j*randn(ib,R);
        Chat = randn(ic,R) + 1j*randn(ic,R);
        aux = 1000;
        error = zeros(1,aux);
        error(1) = ((norm((mode_1 - Ahat*(tensor.mtx_prod_kr(Chat,Bhat).')),'fro'))^2)/((norm(mode_1,'fro')^2));
        for i = 2:aux
            Bhat = mode_2*pinv((tensor.mtx_prod_kr(Chat,Ahat)).');
            Chat = mode_3*pinv((tensor.mtx_prod_kr(Bhat,Ahat)).');
            Ahat = mode_1*pinv((tensor.mtx_prod_kr(Chat,Bhat)).');
            error(i) = ((norm((mode_1 - Ahat*(tensor.mtx_prod_kr(Chat,Bhat).')),'fro'))^2)/((norm(mode_1,'fro')^2));
            if abs(error(i) - error(i-1)) < 1e-6
                error = error(1:i);
                break;
            else
                continue;
            end
        end
        Xhat = tensor.fold(Ahat*(tensor.mtx_prod_kr(Chat,Bhat).'),[I J K],1);
        
        aux1 = (norm(A- Ahat,'fro')^2)/(norm(A,'fro')^2);
        nmse1(snr,1) = nmse1(snr,1) + 20*log10(aux1);
        aux2 = (norm(B- Bhat,'fro')^2)/(norm(B,'fro')^2);
        nmse2(snr,1) = nmse2(snr,1) + 20*log10(aux2);
        aux3 = (norm(C- Chat,'fro')^2)/(norm(C,'fro')^2);
        nmse3(snr,1) = nmse3(snr,1) + 20*log10(aux3);
        aux4 = (norm(tensor.unfold(X- Xhat,1),'fro')^2)/(norm(tensor.unfold(X,1),'fro')^2);
        nmse4(snr,1) = nmse4(snr,1) + 20*log10(aux4);
    end
end
nmse1  = nmse1/1000;
nmse2  = nmse2/1000;
nmse3  = nmse3/1000;
nmse4  = nmse4/1000;

figure
txt = ['\bf A'];
plot(SNR,nmse1,'--d','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['\bf B'];
plot(SNR,nmse2,'--d','color', [0.8500 0.3250 0.0980], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['\bf C'];
plot(SNR,nmse3,'--d','color', [0.4660 0.6740 0.1880], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['Reconstruction'];
plot(SNR,nmse4,'-o','color', [0 0.4470 0.7410], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title(['ALS performance: I = ' num2str(I), ', J = ' num2str(J), ', K = ' num2str(K), ' and R = ' num2str(R)])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "southwest");
set(legend_copy,'Interpreter','tex','location','southwest',"fontsize", 12)
grid on;
saveas(gcf,'hw11a2.png')
