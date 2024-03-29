%% ----- Homework 10 ----- %%
clc;
clear;
close all;

load('homework_10_MLSKronF.mat')

rows = [4 4 6];
columns = [3 2 5];

%% Initialization by HOSVD
[Matrices] = tensor.MLSKronF(X,rows,columns,'1');
aux1 = tensor.mtx_prod_kron(Matrices{1},Matrices{2});
Xhat = tensor.mtx_prod_kron(aux1,Matrices{3});
disp(['Checking the NMSE (dB) between the original matrix X and its'...
    'reconstruction with MLSKronF:'])
nmse1 = (norm(X - Xhat,'fro')^2)/(norm(X,'fro')^2);
nmse1 = 20*log10(nmse1)
disp(['Checking the NMSE (dB) between the original matrix A and its'...
    'estimation:'])
nmsea = (norm(A - Matrices{1},'fro')^2)/(norm(A,'fro')^2);
nmsea = 20*log10(nmsea)
disp(['Checking the NMSE (dB) between the original matrix B and its'...
    'estimation:'])
nmseb = (norm(B - Matrices{2},'fro')^2)/(norm(B,'fro')^2);
nmseb = 20*log10(nmseb)
disp(['Checking the NMSE (dB) between the original matrix C and its'...
    'estimation:'])
nmsec = (norm(C - Matrices{3},'fro')^2)/(norm(C,'fro')^2);
nmsec = 20*log10(nmsec)

%% Initialization by HOOI
[Matrices] = tensor.MLSKronF(X,rows,columns,'2');
aux1 = tensor.mtx_prod_kron(Matrices{1},Matrices{2});
Xhat = tensor.mtx_prod_kron(aux1,Matrices{3});
disp(['Checking the NMSE (dB) between the original matrix X and its'...
    'reconstruction with MLSKronF:'])
nmse1 = (norm(X - Xhat,'fro')^2)/(norm(X,'fro')^2);
nmse1 = 20*log10(nmse1)
disp(['Checking the NMSE (dB) between the original matrix A and its'...
    'estimation:'])
nmsea = (norm(A - Matrices{1},'fro')^2)/(norm(A,'fro')^2);
nmsea = 20*log10(nmsea)
disp(['Checking the NMSE (dB) between the original matrix B and its'...
    'estimation:'])
nmseb = (norm(B - Matrices{2},'fro')^2)/(norm(B,'fro')^2);
nmseb = 20*log10(nmseb)
disp(['Checking the NMSE (dB) between the original matrix C and its'...
    'estimation:'])
nmsec = (norm(C - Matrices{3},'fro')^2)/(norm(C,'fro')^2);
nmsec = 20*log10(nmsec)

%% Monte Carlo Experiment
ia = 2;
ib = 2;
ic = 2;
ja = 2;
jb = 2;
jc = 2;
rows = [ia ib ic];
columns = [ja jb jc];
SNR = [0 5 10 15 20 25 30];
nmse1 = zeros(length(SNR),1);
nmse2 = zeros(length(SNR),1);
for snr = 1:length(SNR)
    for mc = 1:1000
        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise/2)*(randn(ia*ib*ic,ja*jb*jc) ...
            + 1j*randn(ia*ib*ic,ja*jb*jc));
        
        A = randn(ia,ja) + 1i*randn(ia,ja);
        B = randn(ib,jb) + 1i*randn(ib,jb);
        C = randn(ic,jc) + 1i*randn(ic,jc);
        X = tensor.mtx_prod_kron(A,B);
        X = tensor.mtx_prod_kron(X,C);
        X_noisy = X + noise;
        
        % HOSVD
        [Matrices] = tensor.MLSKronF(X_noisy,rows,columns,'1');
        aux1 = tensor.mtx_prod_kron(Matrices{1},Matrices{2});
        Xhat = tensor.mtx_prod_kron(aux1,Matrices{3});
        aux1 = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse1(snr,1) = nmse1(snr,1) + 20*log10(aux1);
        % HOOI
        [Matrices] = tensor.MLSKronF(X_noisy,rows,columns,'2');
        aux2 = tensor.mtx_prod_kron(Matrices{1},Matrices{2});
        Xhat = tensor.mtx_prod_kron(aux2,Matrices{3});
        aux2 = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse2(snr,1) = nmse2(snr,1) + 20*log10(aux2);
    end
end
nmse1  = nmse1/1000;
nmse2  = nmse2/1000;

figure
txt = ['HOSVD Start'];
plot(SNR,nmse1,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2,...
    "markersize", 8, "DisplayName", txt);
hold on;
txt = ['HOOI Start'];
plot(SNR,nmse2,'-x','color', [0 0.4470 0.7410], "linewidth", 2,...
    "markersize", 8, "DisplayName", txt);
hold off;
title(['MLSKronF NMSE vs. SNR: (2,2,2,2,2,2)'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
grid on;
saveas(gcf,'hw10a1.png')

ia = 5;
ib = 5;
ic = 5;
ja = 5;
jb = 5;
jc = 5;
rows = [ia ib ic];
columns = [ja jb jc];
SNR = [0 5 10 15 20 25 30];
nmse1 = zeros(length(SNR),1);
nmse2 = zeros(length(SNR),1);
for snr = 1:length(SNR)
    for mc = 1:250
        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise/2)*(randn(ia*ib*ic,ja*jb*jc) ...
            + 1j*randn(ia*ib*ic,ja*jb*jc));
        
        A = randn(ia,ja) + 1i*randn(ia,ja);
        B = randn(ib,jb) + 1i*randn(ib,jb);
        C = randn(ic,jc) + 1i*randn(ic,jc);
        X = tensor.mtx_prod_kron(A,B);
        X = tensor.mtx_prod_kron(X,C);
        X_noisy = X + noise;
        
        % HOSVD
        [Matrices] = tensor.MLSKronF(X_noisy,rows,columns,'1');
        aux1 = tensor.mtx_prod_kron(Matrices{1},Matrices{2});
        Xhat = tensor.mtx_prod_kron(aux1,Matrices{3});
        aux1 = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse1(snr,1) = nmse1(snr,1) + 20*log10(aux1);
        % HOOI
        [Matrices] = tensor.MLSKronF(X_noisy,rows,columns,'2');
        aux2 = tensor.mtx_prod_kron(Matrices{1},Matrices{2});
        Xhat = tensor.mtx_prod_kron(aux2,Matrices{3});
        aux2 = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse2(snr,1) = nmse2(snr,1) + 20*log10(aux2);
    end
end
nmse1  = nmse1/250;
nmse2  = nmse2/250;

figure
txt = ['HOSVD Start'];
plot(SNR,nmse1,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2,...
    "markersize", 8, "DisplayName", txt);
hold on;
txt = ['HOOI Start'];
plot(SNR,nmse2,'-x','color', [0 0.4470 0.7410], "linewidth", 2,...
    "markersize", 8, "DisplayName", txt);
hold off;
title(['MLSKronF NMSE vs. SNR: (5,5,5,5,5,5)'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
grid on;
saveas(gcf,'hw10a2.png')

ia = 2;
ib = 3;
ic = 4;
ja = 2;
jb = 3;
jc = 4;
rows = [ia ib ic];
columns = [ja jb jc];
SNR = [0 5 10 15 20 25 30];
nmse1 = zeros(length(SNR),1);
nmse2 = zeros(length(SNR),1);
for snr = 1:length(SNR)
    for mc = 1:1000
        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise/2)*(randn(ia*ib*ic,ja*jb*jc) ...
            + 1j*randn(ia*ib*ic,ja*jb*jc));
        
        A = randn(ia,ja) + 1i*randn(ia,ja);
        B = randn(ib,jb) + 1i*randn(ib,jb);
        C = randn(ic,jc) + 1i*randn(ic,jc);
        X = tensor.mtx_prod_kron(A,B);
        X = tensor.mtx_prod_kron(X,C);
        X_noisy = X + noise;
        
        % HOSVD
        [Matrices] = tensor.MLSKronF(X_noisy,rows,columns,'1');
        aux1 = tensor.mtx_prod_kron(Matrices{1},Matrices{2});
        Xhat = tensor.mtx_prod_kron(aux1,Matrices{3});
        aux1 = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse1(snr,1) = nmse1(snr,1) + 20*log10(aux1);
        % HOOI
        [Matrices] = tensor.MLSKronF(X_noisy,rows,columns,'2');
        aux2 = tensor.mtx_prod_kron(Matrices{1},Matrices{2});
        Xhat = tensor.mtx_prod_kron(aux2,Matrices{3});
        aux2 = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse2(snr,1) = nmse2(snr,1) + 20*log10(aux2);
    end
end
nmse1  = nmse1/1000;
nmse2  = nmse2/1000;

figure
txt = ['HOSVD Start'];
plot(SNR,nmse1,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2,...
    "markersize", 8, "DisplayName", txt);
hold on;
txt = ['HOOI Start'];
plot(SNR,nmse2,'-x','color', [0 0.4470 0.7410], "linewidth", 2,...
    "markersize", 8, "DisplayName", txt);
hold off;
title(['MLSKronF NMSE vs. SNR: (2,3,4,2,3,4)'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
grid on;
saveas(gcf,'hw10a3.png')

ia = 2;
ib = 3;
ic = 4;
ja = 5;
jb = 6;
jc = 7;
rows = [ia ib ic];
columns = [ja jb jc];
SNR = [0 5 10 15 20 25 30];
nmse1 = zeros(length(SNR),1);
nmse2 = zeros(length(SNR),1);
for snr = 1:length(SNR)
    for mc = 1:1000
        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise/2)*(randn(ia*ib*ic,ja*jb*jc) ...
            + 1j*randn(ia*ib*ic,ja*jb*jc));
        
        A = randn(ia,ja) + 1i*randn(ia,ja);
        B = randn(ib,jb) + 1i*randn(ib,jb);
        C = randn(ic,jc) + 1i*randn(ic,jc);
        X = tensor.mtx_prod_kron(A,B);
        X = tensor.mtx_prod_kron(X,C);
        X_noisy = X + noise;
        
        % HOSVD
        [Matrices] = tensor.MLSKronF(X_noisy,rows,columns,'1');
        aux1 = tensor.mtx_prod_kron(Matrices{1},Matrices{2});
        Xhat = tensor.mtx_prod_kron(aux1,Matrices{3});
        aux1 = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse1(snr,1) = nmse1(snr,1) + 20*log10(aux1);
        % HOOI
        [Matrices] = tensor.MLSKronF(X_noisy,rows,columns,'2');
        aux2 = tensor.mtx_prod_kron(Matrices{1},Matrices{2});
        Xhat = tensor.mtx_prod_kron(aux2,Matrices{3});
        aux2 = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse2(snr,1) = nmse2(snr,1) + 20*log10(aux2);
    end
end
nmse1  = nmse1/1000;
nmse2  = nmse2/1000;

figure
txt = ['HOSVD Start'];
plot(SNR,nmse1,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2,...
    "markersize", 8, "DisplayName", txt);
hold on;
txt = ['HOOI Start'];
plot(SNR,nmse2,'-x','color', [0 0.4470 0.7410], "linewidth", 2,...
    "markersize", 8, "DisplayName", txt);
hold off;
title(['MLSKronF NMSE vs. SNR: (2,3,4,5,6,7)'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
grid on;
saveas(gcf,'hw10a4.png')