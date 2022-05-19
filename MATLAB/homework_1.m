%% ----- Homework 1 ----- %%
clc;
clear;
close all;

N = [2 4 8 16 32 64 128];
time1 = zeros(length(N),1);
time2 = zeros(length(N),1);
for nn = 1:length(N)
    for mc = 1:1000
        A = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
        B = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
        tic;
        A.*B;
        aux = toc;
        time1(nn,1) = time1(nn,1) + aux;
        tic;
        tensor.mtx_prod_had(A,B);
        aux = toc;
        time2(nn,1) = time2(nn,1) + aux;
    end
end
time1 = time1/1000;
time2 = time2/1000;

figure
txt = ['Native Function'];
semilogy(N,time1,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['Author Function'];
semilogy(N,time2,'-o','color', [0.8500 0.3250 0.0980], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title(['Hadamard product performance'])
xlabel('Number of columns')
ylabel('Time (s)')
legend_copy = legend("location", "southeast");
set(legend_copy,'Interpreter','tex','location','southeast',"fontsize", 12)
grid on;
saveas(gcf,'hw1a1.png')

N = [2 4 8 16 32 64];
time1 = zeros(length(N),1);
time2 = zeros(length(N),1);
for nn = 1:length(N)
    for mc = 1:1000
        A = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
        B = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
        tic;
        kron(A,B);
        aux = toc;
        time1(nn,1) = time1(nn,1) + aux;
        tic;
        tensor.mtx_prod_kron(A,B);
        aux = toc;
        time2(nn,1) = time2(nn,1) + aux;
    end
end
time1 = time1/1000;
time2 = time2/1000;

figure
txt = ['Native Function'];
semilogy(N,time1,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['Author Function'];
semilogy(N,time2,'-o','color', [0.8500 0.3250 0.0980], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title(['Kronecker product performance'])
xlabel('Number of columns')
ylabel('Time (s)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
grid on;
saveas(gcf,'hw1a2.png')

N = [2 4 8 16 32 64];
time1 = zeros(length(N),1);
time2 = zeros(length(N),1);
for nn = 1:length(N)
    for mc = 1:1000
        A = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
        B = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
        tic;
        kr(A,B);
        aux = toc;
        time1(nn,1) = time1(nn,1) + aux;
        tic;
        tensor.mtx_prod_kr(A,B);
        aux = toc;
        time2(nn,1) = time2(nn,1) + aux;
    end
end
time1 = time1/1000;
time2 = time2/1000;

figure
txt = ['Native Function'];
semilogy(N,time1,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['Author Function'];
semilogy(N,time2,'-o','color', [0.8500 0.3250 0.0980], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title(['Khatri-Rao product performance'])
xlabel('Number of columns')
ylabel('Time (s)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
grid on;
saveas(gcf,'hw1a3.png')
