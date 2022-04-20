%% ----- Homework 9 ----- %%
clc;
clear;
close all;
I = 3;
J = 5;
K = 7;
R = 3;

A = randn(I,R) + 1j*randn(I,R);
B = randn(J,R) + 1j*randn(J,R);
C = randn(K,R) + 1j*randn(K,R);
X = tensor.fold(A*(tensor.mtx_prod_kr(C,B).'),[I J K],1);

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


figure
txt = ['\bf Convergence between iterations'];
plot(1:i,20*log10(error),'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
title(['ALS convergence: I = ' num2str(I), ', J = ' num2str(J), ', K = ' num2str(K), ' and R = ' num2str(R)])
xlabel('Iterations')
ylabel('Error (dB)')
legend_copy = legend("location", "southwest");
set(legend_copy,'Interpreter','tex','location','southwest',"fontsize", 12)
grid on;
saveas(gcf,'hw9a1.png')

I = 10;
J = 4;
K = 2;
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
plot(SNR,nmse1,'-d','color', [1 0 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['\bf B'];
plot(SNR,nmse2,'-o','color', [0 1 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['\bf C'];
plot(SNR,nmse3,'-x','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['Reconstruction'];
plot(SNR,nmse4,'-p','color', [0 0 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title(['ALS performance: I = ' num2str(I), ', J = ' num2str(J), ', K = ' num2str(K), ' and R = ' num2str(R)])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "southwest");
set(legend_copy,'Interpreter','tex','location','southwest',"fontsize", 12)
grid on;
saveas(gcf,'hw9a2.png')