%% ----- Homework 13 ----- %%
clc;
clear;
close all;

%% Tensor Train SVD for a fourth order tensor
I1 = 5;
I2 = 5;
I3 = 5;
I4 = 5;
Ranks = [3 3 3];

G1 = randn(I1,Ranks(1));
G2 = randn(Ranks(1),I2,Ranks(2));
G3 = randn(Ranks(2),I3,Ranks(3));
G4 = randn(Ranks(3),I4);

tenX = tensor.contraction(G1,2,G2,1);
tenX = tensor.contraction(tenX,3,G3,1);
tenX = tensor.contraction(tenX,4,G4,1);

[G] = tensor.TTSVD(tenX,Ranks);

tenXhat = tensor.contraction(G{1},2,G{2},1);
tenXhat = tensor.contraction(tenXhat,3,G{3},1);
tenXhat = tensor.contraction(tenXhat,4,G{4},1);

disp('Checking the TTSVD NMSE (dB) from reconstruction using ranks = [3 3 3]:')
nmsex = (norm(tensor.unfold(tenX- tenXhat,1),'fro')^2)/(norm(tensor.unfold(tenX,1),'fro')^2);
nmsex = 20*log10(nmsex)

Ranks = [2 2 2];
[G] = tensor.TTSVD(tenX,Ranks);

tenXhat = tensor.contraction(G{1},2,G{2},1);
tenXhat = tensor.contraction(tenXhat,3,G{3},1);
tenXhat = tensor.contraction(tenXhat,4,G{4},1);
disp('Checking the TTSVD NMSE (dB) from reconstruction using ranks = [2 2 2]:')
nmsex = (norm(tensor.unfold(tenX- tenXhat,1),'fro')^2)/(norm(tensor.unfold(tenX,1),'fro')^2);
nmsex = 20*log10(nmsex)

%% Monte Carlo Simulation
I1 = 5;
I2 = 5;
I3 = 5;
I4 = 5;
Ranks  = [3 3 3];
Ranks1 = [2 2 2];
Ranks2 = [3 3 3];

SNR = [0 5 10 15 20 25 30];
nmse1 = zeros(length(SNR),1);
nmse2 = zeros(length(SNR),1);
for snr = 1:length(SNR)
    snr
    for mc = 1:1000
        G1 = randn(I1,Ranks(1));
        G2 = randn(Ranks(1),I2,Ranks(2));
        G3 = randn(Ranks(2),I3,Ranks(3));
        G4 = randn(Ranks(3),I4);

        tenX = tensor.contraction(G1,2,G2,1);
        tenX = tensor.contraction(tenX,3,G3,1);
        tenX = tensor.contraction(tenX,4,G4,1);

        var_noise = 1/(10^(SNR(snr)/10));
        noise = sqrt(var_noise)*(randn(size(tenX)));
        tenX_noisy = tenX + noise;
        
        [GR1] = tensor.TTSVD(tenX_noisy,Ranks1);
        tenXhat = tensor.contraction(GR1{1},2,GR1{2},1);
        tenXhat = tensor.contraction(tenXhat,3,GR1{3},1);
        tenXhat = tensor.contraction(tenXhat,4,GR1{4},1);
        nmse1(snr,1) = nmse1(snr,1) + 20*log10((norm(tensor.unfold(tenX- tenXhat,1),'fro')^2)/(norm(tensor.unfold(tenX,1),'fro')^2));
        
        [GR2] = tensor.TTSVD(tenX_noisy,Ranks2);
        tenXhat = tensor.contraction(GR2{1},2,GR2{2},1);
        tenXhat = tensor.contraction(tenXhat,3,GR2{3},1);
        tenXhat = tensor.contraction(tenXhat,4,GR2{4},1);
        nmse2(snr,1) = nmse2(snr,1) + 20*log10((norm(tensor.unfold(tenX- tenXhat,1),'fro')^2)/(norm(tensor.unfold(tenX,1),'fro')^2));
        
    end
end
nmse1  = nmse1/1000;
nmse2  = nmse2/1000;

figure
txt = ['Ranks = [2 2 2]'];
plot(SNR,nmse1,'--d','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['Ranks = [3 3 3]'];
plot(SNR,nmse2,'-o','color', [0 0.4470 0.7410], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title(['TTSVD performance under imperfect scenario'])
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
legend_copy = legend("location", "southwest");
set(legend_copy,'Interpreter','tex',"fontsize", 12)
grid on;
saveas(gcf,'hw13.png')