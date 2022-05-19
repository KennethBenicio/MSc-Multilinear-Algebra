%% ----- Homework 0 ----- %%
clc;
clear;
close all;

N = [1 2 3 4 5 6 7 8];
time1 = zeros(length(N),1);
time2 = zeros(length(N),1);
for nn = 1:length(N)
    for mc = 1:5000
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
time1 = time1/5000;
time2 = time2/5000;

figure
txt = ['(\bf A_{N \times N} \otimes B_{N \times N})^{-1}'];
semilogy(N,time1,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['\bf A^{-1}_{N \times N} \otimes B^{-1}_{N \times N}'];
semilogy(N,time2,'-o','color', [0.8500 0.3250 0.0980], "linewidth", 2, "markersize", 8, "DisplayName", txt);
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
    for mc = 1:1000
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
                A1 = tensor.mtx_prod_kron(inv(A1),inv(A2));
            end    
        end
        aux = toc;
        time2(kk,1) = time2(kk,1) + aux;
    end
end
time1 = time1/1000;
time2 = time2/1000;

figure
txt = ['\bf (\otimes^{K}_{i = 1}A^{(i)}_{4 \times 4})^{-1}'];
semilogy(K,time1,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['\bf \otimes^{K}_{i = 1}(A^{(i)}_{4 \times 4})^{-1}'];
semilogy(K,time2,'-o','color', [0.8500 0.3250 0.0980], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title(['Multiple Kronecker products performance'])
xlabel('Number of products')
ylabel('Time (s)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
grid on;
saveas(gcf,'hw0a2.png')