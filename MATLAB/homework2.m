%% ----- Homework 2 ----- %%
clc;
clear;
close all;

R = 2;
I = [2 4 8 16 32 64];
time1 = zeros(length(I),1);
time2 = zeros(length(I),1);
time3 = zeros(length(I),1);
for ii = 1:length(I)
    for mc = 1:250
        A = randn(I(ii),I(ii)) + 1j*randn(I(ii),I(ii));
        B = randn(I(ii),I(ii)) + 1j*randn(I(ii),I(ii));
        tic;
        pinv(tensor.mtx_prod_kr(A,B));
        aux = toc;
        time1(ii,1) = time1(ii,1) + aux;
        tic;
        (tensor.mtx_prod_kr(A,B).'*tensor.mtx_prod_kr(A,B))\(tensor.mtx_prod_kr(A,B).');
        aux = toc;
        time2(ii,1) = time2(ii,1) + aux;
        tic;
        tensor.mtx_prod_had((A.'*A),(B.'*B))\(tensor.mtx_prod_kr(A,B).');
        aux = toc;
        time3(ii,1) = time3(ii,1) + aux;
    end
end
time1 = time1/250;
time2 = time2/250;
time3 = time3/250;

figure
txt = ['1st Method'];
plot(I,time1,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['2nd Method'];
plot(I,time2,'-d','color', [0.8500 0.3250 0.0980], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['3rd Method'];
plot(I,time3,'-d','color', [0.4660 0.6740 0.1880], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title(['Khatri-Rao product performance considering 2 columns'])
xlabel('Number of rows')
ylabel('Time (s)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
grid on;
saveas(gcf,'hw2a1.png')

R = 4;
I = [2 4 8 16 32 64];
time1 = zeros(length(I),1);
time2 = zeros(length(I),1);
time3 = zeros(length(I),1);
for ii = 1:length(I)
    for mc = 1:250
        A = randn(I(ii),I(ii)) + 1j*randn(I(ii),I(ii));
        B = randn(I(ii),I(ii)) + 1j*randn(I(ii),I(ii));
        tic;
        pinv(tensor.mtx_prod_kr(A,B));
        aux = toc;
        time1(ii,1) = time1(ii,1) + aux;
        tic;
        (tensor.mtx_prod_kr(A,B).'*tensor.mtx_prod_kr(A,B))\(tensor.mtx_prod_kr(A,B).');
        aux = toc;
        time2(ii,1) = time2(ii,1) + aux;
        tic;
        tensor.mtx_prod_had((A.'*A),(B.'*B))\(tensor.mtx_prod_kr(A,B).');
        aux = toc;
        time3(ii,1) = time3(ii,1) + aux;
    end
end
time1 = time1/250;
time2 = time2/250;
time3 = time3/250;

figure
txt = ['1st Method'];
plot(I,time1,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['2nd Method'];
plot(I,time2,'-d','color', [0.8500 0.3250 0.0980], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['3rd Method'];
plot(I,time3,'-d','color', [0.4660 0.6740 0.1880], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title(['Khatri-Rao product performance considering 4 columns'])
xlabel('Number of rows')
ylabel('Time (s)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
grid on;
saveas(gcf,'hw2a2.png')

I = 4;
R = 2;
N = [2 4 6 8];
time1 = zeros(length(N),1);
time2 = zeros(length(N),1);
for nn = 1:length(N)
    for mc = 1:250 
        tic;
        for ii = 1:N(nn)
            if ii == 1
                A1 = randn(I,R) + 1j*randn(I,R); 
                continue
            else
                A2 = randn(I,R) + 1j*randn(I,R);
                A1 = tensor.mtx_prod_kron(A1,A2);
            end    
        end
        aux = toc;
        time1(nn,1) = time1(nn,1) + aux;
    end
end
time1 = time1/250;

figure
txt = ['\bf ⋄^{N}_{n = 1} (A^{(n)}_{4 \times 2})'];
semilogy(N,time1,'-d','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
title(['Multiple Khatri-Rao products performance'])
xlabel('Number of Products')
ylabel('Time (s)')
legend_copy = legend("location", "northwest");
set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
grid on;
saveas(gcf,'hw2a3.png')
