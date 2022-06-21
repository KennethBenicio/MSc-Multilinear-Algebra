%% Tensor LMS Algorithm
clc;
close all;
clear all;

% System parameters
Nh = 5;
Nv = 5;
QAM = 2;
runs = 1000;
buffer = Nh*Nv;
Samples = 15000;

% Noise energy
snr = 30;
var_noise = 1/(10^(snr/10));  

%% 0.1 Learning Step
mu = 0.1;

% Monter Carlo simulations
mse = zeros(Samples,1);
for rr = 1:runs
    % Filter weights
    weights_v = complex(randn(Nv,Samples),randn(Nv,Samples));
    weights_h = complex(randn(Nh,Samples),randn(Nh,Samples));
    
    % M-QAM source signal
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
   
    % Adding noise to the M-QAM source signal
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    
    % TLMS algorithm
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
txt = ['\mu = 0.1'];
plot(1:Samples, 10*log10(mse),'-','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;

%% 0.5 Learning Step
mu = 0.5;

mse = zeros(Samples,1);
for rr = 1:runs
    % Filter weights
    weights_v = complex(randn(Nv,Samples),randn(Nv,Samples));
    weights_h = complex(randn(Nh,Samples),randn(Nh,Samples));
    
    % M-QAM source signal
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
   
    % Adding noise to the M-QAM source signal
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    
    % TLMS algorithm
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

txt = ['\mu = 0.5'];
plot(1:Samples, 10*log10(mse),'-','color', [0.8500 0.3250 0.0980], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;

%% 0.05 Learning Step
mu = 0.05;

mse = zeros(Samples,1);
for rr = 1:runs
    % Filter weights
    weights_v = complex(randn(Nv,Samples),randn(Nv,Samples));
    weights_h = complex(randn(Nh,Samples),randn(Nh,Samples));
    
    % M-QAM source signal
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
   
    % Adding noise to the M-QAM source signal
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    
    % TLMS algorithm
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

txt = ['\mu = 0.05'];
plot(1:Samples, 10*log10(mse),'-','color', [0.4660 0.6740 0.1880], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;

%% 0.005 Learning Step
mu = 0.005;

mse = zeros(Samples,1);
for rr = 1:runs
    % Filter weights
    weights_v = complex(randn(Nv,Samples),randn(Nv,Samples));
    weights_h = complex(randn(Nh,Samples),randn(Nh,Samples));
    
    % M-QAM source signal
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
   
    % Adding noise to the M-QAM source signal
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    
    % TLMS algorithm
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

txt = ['\mu = 0.005'];
plot(1:Samples, 10*log10(mse),'-','color', [0.9290 0.6940 0.1250], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
title('TLMS Behavior for Different Learning Steps');
xlabel('Iteration');
ylabel('MSE');
legend_copy = legend("location", "northeast");
set (legend_copy, "fontsize", 12);
grid on;
saveas(gcf,'tlms_mse.png')