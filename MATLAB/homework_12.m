%% ----- Homework 12 ----- %%
clc;
clear;
close all;

%% TKPSVD 
R = 1;
N = 3;

A = zeros(80,R);
B = zeros(20,R);
C = zeros(8,R);
tenX = zeros(100,16,8);
for r = 1:R
    tenA = randn(10,4,2);
    %tenA = tenA./frob(tenA);
    A(:,r)= tenA(:);
    tenB = randn(5,2,2);
    %tenB = tenB./frob(tenB);
    B(:,r)= tenB(:);
    tenC = randn(2,2,2);
    %tenC = tenC./frob(tenC);
    C(:,r)= tenC(:);
    varr = tensor.ten_prod_kron(tenC,tenB);
    tenX = tenX + tensor.ten_prod_kron(varr,tenA);
end

aux{1} = num2cell(size(tenA));
aux{2} = num2cell(size(tenB));
aux{3} = num2cell(size(tenC));
tenSize = horzcat(aux{:});
tenDim = cellfun(@(aux) prod(cell2mat(aux)), aux,'UniformOutput',false);

tenX = tensor.ten_prod_kron(tenC,tenB);
tenX = tensor.ten_prod_kron(tenX,tenA);

[Ahat,Bhat,Chat] = tensor.TKPSVD(tenX,tenSize,tenDim,N,R);

aux1 = aux{1};
aux2 = aux{2};
aux3 = aux{3};
tenXhat = zeros(100,16,8);
for r = 1:R
    tenAhat_r = reshape(Ahat(:,r),aux1{:});
    %tenAhat_r = tenAhat_r./frob(tenAhat_r);
    tenBhat_r = reshape(Bhat(:,r),aux2{:});
    %tenBhat_r = tenBhat_r./frob(tenBhat_r);
    tenChat_r = reshape(Chat(:,r),aux3{:});
    %tenChat_r = tenChat_r./frob(tenChat_r);
    varr = tensor.ten_prod_kron(tenChat_r,tenBhat_r);
    tenXhat = tenXhat + tensor.ten_prod_kron(varr,tenAhat_r);
end

disp('Checking the NMSE (dB) between the original tensor X and its'...
    'reconstruction with TKPSVD:')
nmsex = (norm(tensor.unfold(tenX- tenXhat,1),'fro')^2)...
    /(norm(tensor.unfold(tenX,1),'fro')^2);
nmsex = 20*log10(nmsex)
disp('Checking the NMSE (dB) between the original matrix A and its'...
    'reconstruction with TKPSVD:')
nmsea = 20*log10((norm(A- Ahat,'fro')^2)/(norm(A,'fro')^2))
disp('Checking the NMSE (dB) between the original matrix B and its'...
    'reconstruction with TKPSVD:')
nmseb = 20*log10((norm(B- Bhat,'fro')^2)/(norm(B,'fro')^2))
disp('Checking the NMSE (dB) between the original matrix C and its'...
    'reconstruction with TKPSVD:')
nmsec = 20*log10((norm(C- Chat,'fro')^2)/(norm(C,'fro')^2))

%% Monte Carlo Simulation

R = 2;
N = 3;
SNR = [0 5 10 15 20 25 30];
nmse1 = zeros(length(SNR),1);
nmse2 = zeros(length(SNR),1);
nmse3 = zeros(length(SNR),1);
nmse4 = zeros(length(SNR),1);
for snr = 1:length(SNR)
    snr
    for mc = 1:1000
        A = zeros(80,R);
        B = zeros(20,R);
        C = zeros(8,R);
        tenX = zeros(100,16,8);
        for r = 1:R
            tenA = randn(10,4,2);
            %tenA = tenA./frob(tenA);
            A(:,r)= tenA(:);
            tenB = randn(5,2,2);
            %tenB = tenB./frob(tenB);
            B(:,r)= tenB(:);
            tenC = randn(2,2,2);
            %tenC = tenC./frob(tenC);
            C(:,r)= tenC(:);
            varr = tensor.ten_prod_kron(tenC,tenB);
            tenX = tenX + tensor.ten_prod_kron(varr,tenA);
        end

        aux{1} = num2cell(size(tenA));
        aux{2} = num2cell(size(tenB));
        aux{3} = num2cell(size(tenC));
        tenSize = horzcat(aux{:});
        tenDim = cellfun(@(aux) prod(cell2mat(aux)), aux,...
            'UniformOutput',false);

        tenX = tensor.ten_prod_kron(tenC,tenB);
        tenX = tensor.ten_prod_kron(tenX,tenA);
        
        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise)*(randn(size(tenX)));
        tenX_noisy = tenX + noise;
        
        [Ahat,Bhat,Chat] = tensor.TKPSVD(tenX_noisy,tenSize,tenDim,N,R);

        aux1 = aux{1};
        aux2 = aux{2};
        aux3 = aux{3};
        tenXhat = zeros(100,16,8);
        for r = 1:R
            tenAhat_r = reshape(Ahat(:,r),aux1{:});
            %tenAhat_r = tenAhat_r./frob(tenAhat_r);
            tenBhat_r = reshape(Bhat(:,r),aux2{:});
            %tenBhat_r = tenBhat_r./frob(tenBhat_r);
            tenChat_r = reshape(Chat(:,r),aux3{:});
            %tenChat_r = tenChat_r./frob(tenChat_r);
            varr = tensor.ten_prod_kron(tenChat_r,tenBhat_r);
            tenXhat = tenXhat + tensor.ten_prod_kron(varr,tenAhat_r);
        end
        
        var1 = (norm(A - Ahat,'fro')^2)/(norm(A,'fro')^2);
        nmse1(snr,1) = nmse1(snr,1) + 20*log10(var1);
        var2 = (norm(B - Bhat,'fro')^2)/(norm(B,'fro')^2);
        nmse2(snr,1) = nmse2(snr,1) + 20*log10(var2);
        var3 = (norm(C - Chat,'fro')^2)/(norm(C,'fro')^2);
        nmse3(snr,1) = nmse3(snr,1) + 20*log10(var3);
        nmse4(snr,1) = nmse4(snr,1) + 20*log10((norm(tensor.unfold(tenX...
            - tenXhat,1),'fro')^2)/(norm(tensor.unfold(tenX,1),'fro')^2));
        
    end
end
nmse1  = nmse1/1000;
nmse2  = nmse2/1000;
nmse3  = nmse3/1000;
nmse4  = nmse4/1000;

figure
txt = ['\bf A'];
plot(SNR,nmse1,'--d','color', [0.3010 0.7450 0.9330], "linewidth", 2,...
    "markersize", 8, "DisplayName", txt);
hold on;
txt = ['\bf B'];
plot(SNR,nmse2,'--d','color', [0.8500 0.3250 0.0980], "linewidth", 2,...
    "markersize", 8, "DisplayName", txt);
hold on;
txt = ['\bf C'];
plot(SNR,nmse3,'--d','color', [0.4660 0.6740 0.1880], "linewidth", 2,...
    "markersize", 8, "DisplayName", txt);
hold on;
txt = ['Reconstruction'];
plot(SNR,nmse4,'-o','color', [0 0.4470 0.7410], "linewidth", 2,...
    "markersize", 8, "DisplayName", txt);
hold off;
title(['TKPSVD performance under imperfect scenario'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "southwest");
set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
grid on;
saveas(gcf,'hw12.png')