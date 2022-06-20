%% Alternating Tensor LMS Algorithm
clc;
close all;
clear all;

Nh = 5;
Nv = 5;
Kh = 10;
Kv = 10;

runs = 5;
mu = 0.005;

QAM = 2;
snr = 30;
var_noise = 1/(10^(snr/10));  

buffer = Nh*Nv;
Samples = 10000;
mse = zeros(Samples,1);
Kb = floor(Samples/(Kh + Kv));
for rr = 1:runs
    
    weights_v = complex(randn(Nv,Samples),randn(Nv,Samples));
    weights_h = complex(zeros(Nh,Samples),randn(Nh,Samples));
    
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
   
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    kh = 1;
    kv = 1;
    for kk = 1:(Kh + Kv):(Kb*(Kh + Kv)) - 2*(Kh + Kv)
        for kh = kk:(kk+Kh-1)
            xmtx = signal_x(kh:kh+buffer-1);
            uh = reshape(xmtx,[Nh Nv])*conj(weights_v(:,kh));
            error = signal_d(kh) - (tensor.mtx_prod_kron(weights_v(:,kh),weights_h(:,kh))')*xmtx;
            mu_norm = mu/(norm(uh,2)^2);
            weights_h(:,kh+1) = weights_h(:,kh) + mu_norm*uh*conj(error); 
            weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1));
            mse(kh,1) =  mse(kh,1) + abs(signal_d(kh) - weights'*xmtx).^2;
        end
        for kv = (kk+Kh):(kk+Kh+Kv-1)
            xmtx = signal_x(kv+Kv:kv+buffer+Kv-1);
            uv = (reshape(xmtx,[Nh Nv]).')*conj(weights_h(:,kv));
            error = signal_d(kv) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1))')*xmtx;
            mu_norm = mu/(norm(uv,2)^2);
            weights_v(:,kv+1) = weights_v(:,kv) + mu_norm*uv*conj(error); 
            weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1));
            mse(kv,1) =  mse(kv,1) + abs(signal_d(kv) - weights'*xmtx).^2;
        end
    end
end
mse = mse/runs;

% MSE Curve
figure
plot(1:Samples, 10*log10(mse),'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 8);
title('ATLMS Behavior');
xlabel('Samples');
ylabel('MSE');
grid on;
saveas(gcf,'atlms_mse.png')