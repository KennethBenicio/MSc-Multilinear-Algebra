%% Tensor LMS Algorithm
clc;
close all;
clear all;

Nh = 5;
Nv = 5;

runs = 5;
mu = 0.5;

QAM = 2;
snr = 30;
var_noise = 1/(10^(snr/10));  

buffer = Nh*Nv;
Samples = 10000;
mse = zeros(Samples,1);
for rr = 1:runs
    
    weights_v = complex(randn(Nv,Samples),randn(Nv,Samples));
    weights_h = complex(zeros(Nh,Samples),randn(Nh,Samples));
    
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
   
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    
    for kk = 1:(Samples - buffer - 1)
        xmtx = signal_x(kk:kk+buffer-1);
        uh = reshape(xmtx,[Nh Nv])*conj(weights_v(:,kk));
        uv = (reshape(xmtx,[Nh Nv]).')*conj(weights_h(:,kk));
        error = signal_d(kk) - (tensor.mtx_prod_kron(weights_v(:,kk),weights_h(:,kk))')*xmtx;
        mu_norm = mu/(norm(uh,2)^2 + norm(uv,2)^2);
        weights_h(:,kk+1) = weights_h(:,kk) + mu_norm*uh*conj(error); 
        weights_v(:,kk+1) = weights_v(:,kk) + mu_norm*uv*conj(error); 
        weights = tensor.mtx_prod_kron(weights_v(:,kk),weights_h(:,kk));
        mse(kk,1) =  mse(kk,1) + abs(signal_d(kk) - weights'*xmtx).^2;
    end
end
mse = mse/runs;

% MSE Curve
figure
plot(1:Samples, 10*log10(mse),'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 8);
title('TLMS Behavior');
xlabel('Samples');
ylabel('MSE');
grid on;
saveas(gcf,'tlms_mse.png')