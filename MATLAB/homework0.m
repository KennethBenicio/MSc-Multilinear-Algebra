%% ----- Homework 0 ----- %%
clc;
clear;
close all;

N = [1 2 3 4 5 6 7 8];
time1 = zeros(length(N),1);
time2 = zeros(length(N),1);
for nn = 1:length(N)
    for mc = 1:1000
        A = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
        B = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
        tic;
        inv(tensor.mtx_prod_kron(A,B));
        aux = toc;
        time1(nn,1) = time1(nn,1) + aux;
        tic;
        tensor.mtx_prod_kron(inv(A),inv(B));
        aux = toc;
        time2(nn,1) = time2(nn,1) + aux;
    end
end
time1 = time1/1000;
time2 = time2/1000;

figure
txt = ['(\bf A_{N \times N} \otimes B_{N \times N})^{-1}'];
semilogy(N,time1,'-d','color', [0 1 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['\bf A^{-1}_{N \times N} \otimes B^{-1}_{N \times N}'];
semilogy(N,time2,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title(['Performance of inverse operation'])
xlabel('Number of columns')
ylabel('Time (s)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
grid on;
saveas(gcf,'hw0a1.png')

N = 2;
K = [1 2 3 4 5 6 7 8 9 10];
time1 = zeros(length(K),1);
time2 = zeros(length(K),1);
for kk = 1:length(K)
    for mc = 1:500
        tic;
        for ii = 1:K(kk)    
            if ii == 1
                A1 = randn(N,N) + 1j*randn(N,N); 
                continue
            else
                A2 = randn(N,N) + 1j*randn(N,N);
                A1 = tensor.mtx_prod_kron(A1,A2);
            end    
        end
        inv(A1);
        aux = toc;
        time1(kk,1) = time1(kk,1) + aux;
        
        tic;
        for ii = 1:K(kk)    
            if ii == 1
                A1 = randn(N,N) + 1j*randn(N,N); 
                continue
            else
                A2 = randn(N,N) + 1j*randn(N,N);
                A1 = tensor.mtx_prod_kron(inv(A1),inv(A2));clc;
clear;
close all;

A = randn(4,2) + 1j*randn(4,2);
B = randn(6,3) + 1j*randn(6,3);
X = tensor.mtx_prod_kron(A,B);
[Ahat,Bhat] = tensor.LSKronF(X,4,2,6,3);
Xhat = tensor.mtx_prod_kron(Ahat,Bhat);

nmsex = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
nmsex = 20*log10(nmsex)
nmsea = (norm(A- Ahat,'fro')^2)/(norm(A,'fro')^2);
nmsea = 20*log10(nmsea)
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
        X = X + noise;
        
        [Ahat,Bhat] = tensor.LSKronF(X,I,P,J,Q);
        Xhat = tensor.mtx_prod_kron(Ahat,Bhat);
        aux = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
        nmse(snr,1) = nmse(snr,1) + 20*log10(aux);
    end
end
nmse  = nmse/1000;

figure
txt = ['I = ' num2str(I), ', J = ' num2str(J), ', P = ' num2str(P), ' and Q = ' num2str(Q)];
plot(SNR,nmse,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
title(['LSKronF performance under imperfect scenario'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
grid on;
saveas(gcf,'hw4a1.png')

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

figure
txt = ['I = ' num2str(I), ', J = ' num2str(J), ', P = ' num2str(P), ' and Q = ' num2str(Q)];
plot(SNR,nmse,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
title(['LSKronF performance under imperfect scenario'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
grid on;
saveas(gcf,'hw4a2.png')

            end    
        end
        aux = toc;
        time2(kk,1) = time2(kk,1) + aux;
    end
end
time1 = time1/500;
time2 = time2/500;

figure
txt = ['\bf (\otimes^{K}_{i = 1}A^{(i)}_{4 \times 4})^{-1}'];
semilogy(K,time1,'-d','color', [0 1 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['\bf \otimes^{K}_{i = 1}(A^{(i)}_{4 \times 4})^{-1}'];
semilogy(K,time2,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title(['Multiple Kronecker products performance'])
xlabel('Number of products')
ylabel('Time (s)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
grid on;
saveas(gcf,'hw0a2.png')