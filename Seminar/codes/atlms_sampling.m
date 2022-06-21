%% Alternating Tensor LMS Algorithm Sampling and Time 
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

%% ATLMS with learning step 0.5 and Kh = Kv = 50
mu = 0.5;
Kh = 50;
Kv = 50;

% Monter Carlo simulations
mse = zeros(Samples,1);
Kb = floor(Samples/(Kh + Kv));
tic();
for rr = 1:runs
    % Filter weights
    weights_v = complex(ones(Nv,Samples),ones(Nv,Samples));
    weights_h = complex(ones(Nh,Samples),ones(Nh,Samples));
    
    % M-QAM source signal
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
    % Adding noise to the M-QAM source signal
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    kh = 1;
    kv = 1;
    % ATLMS algorithm
    for kk = 1:(Kh + Kv):(Kb*(Kh + Kv)) - 5*(Kh + Kv)
        for kh = kk:(kk+Kh-1)
            weights_v(:,kh+1) = weights_v(:,kh);
            xmtx = signal_x(kh:kh+buffer-1);
            uh = reshape(xmtx,[Nh Nv])*conj(weights_v(:,kv));
            error = signal_d(kh) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh))')*xmtx;
            mu_norm = mu/(norm(uh,2)^2);
            weights_h(:,kh+1) = weights_h(:,kh) + mu_norm*uh*conj(error); 
            weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh));
            mse(kh,1) =  mse(kh,1) + abs(signal_d(kh) - weights'*xmtx).^2;
        end
        for kv = (kk+Kh):(kk+Kh+Kv-1)
            weights_h(:,kv+1) = weights_h(:,kv);
            xmtx = signal_x(kv:kv+buffer-1);
            uv = (reshape(xmtx,[Nh Nv]).')*conj(weights_h(:,kh+1));
            error = signal_d(kv) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1))')*xmtx;
            mu_norm = mu/(norm(uv,2)^2);
            weights_v(:,kv+1) = weights_v(:,kv) + mu_norm*uv*conj(error); 
            weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1));
            mse(kv,1) =  mse(kv,1) + abs(signal_d(kv) - weights'*xmtx).^2;
        end
    end
end
aux1 = toc();
aux1 = aux1/runs;
mse = mse/runs;

% MSE Curve
figure
txt = ['\mu = 0.5 with K_h = K_v = 50'];
plot(1:Samples, 10*log10(mse),'-','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;

%% ATLMS with learning step 0.5 and Kh = Kv = 250
mu = 0.5;
Kh = 250;
Kv = 250;

% Monter Carlo simulations
mse = zeros(Samples,1);
Kb = floor(Samples/(Kh + Kv));
tic();
for rr = 1:runs
    % Filter weights
    weights_v = complex(ones(Nv,Samples),ones(Nv,Samples));
    weights_h = complex(ones(Nh,Samples),ones(Nh,Samples));
    
    % M-QAM source signal
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
    % Adding noise to the M-QAM source signal
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    kh = 1;
    kv = 1;
    % ATLMS algorithm
    for kk = 1:(Kh + Kv):(Kb*(Kh + Kv)) - 5*(Kh + Kv)
        for kh = kk:(kk+Kh-1)
            weights_v(:,kh+1) = weights_v(:,kh);
            xmtx = signal_x(kh:kh+buffer-1);
            uh = reshape(xmtx,[Nh Nv])*conj(weights_v(:,kv));
            error = signal_d(kh) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh))')*xmtx;
            mu_norm = mu/(norm(uh,2)^2);
            weights_h(:,kh+1) = weights_h(:,kh) + mu_norm*uh*conj(error); 
            weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh));
            mse(kh,1) =  mse(kh,1) + abs(signal_d(kh) - weights'*xmtx).^2;
        end
        for kv = (kk+Kh):(kk+Kh+Kv-1)
            weights_h(:,kv+1) = weights_h(:,kv);
            xmtx = signal_x(kv:kv+buffer-1);
            uv = (reshape(xmtx,[Nh Nv]).')*conj(weights_h(:,kh+1));
            error = signal_d(kv) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1))')*xmtx;
            mu_norm = mu/(norm(uv,2)^2);
            weights_v(:,kv+1) = weights_v(:,kv) + mu_norm*uv*conj(error); 
            weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1));
            mse(kv,1) =  mse(kv,1) + abs(signal_d(kv) - weights'*xmtx).^2;
        end
    end
end
aux2 = toc();
aux2 = aux2/runs;
mse = mse/runs;

txt = ['\mu = 0.5 with K_h = K_v = 250'];
plot(1:Samples, 10*log10(mse),'-','color', [0.4660 0.6740 0.1880], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;

%% ATLMS with learning step 0.5 and Kh = Kv = 500
mu = 0.5;
Kh = 500;
Kv = 500;

% Monter Carlo simulations
mse = zeros(Samples,1);
Kb = floor(Samples/(Kh + Kv));
tic();
for rr = 1:runs
    % Filter weights
    weights_v = complex(ones(Nv,Samples),ones(Nv,Samples));
    weights_h = complex(ones(Nh,Samples),ones(Nh,Samples));
    
    % M-QAM source signal
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
    % Adding noise to the M-QAM source signal
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    kh = 1;
    kv = 1;
    % ATLMS algorithm
    for kk = 1:(Kh + Kv):(Kb*(Kh + Kv)) - 5*(Kh + Kv)
        for kh = kk:(kk+Kh-1)
            weights_v(:,kh+1) = weights_v(:,kh);
            xmtx = signal_x(kh:kh+buffer-1);
            uh = reshape(xmtx,[Nh Nv])*conj(weights_v(:,kv));
            error = signal_d(kh) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh))')*xmtx;
            mu_norm = mu/(norm(uh,2)^2);
            weights_h(:,kh+1) = weights_h(:,kh) + mu_norm*uh*conj(error); 
            weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh));
            mse(kh,1) =  mse(kh,1) + abs(signal_d(kh) - weights'*xmtx).^2;
        end
        for kv = (kk+Kh):(kk+Kh+Kv-1)
            weights_h(:,kv+1) = weights_h(:,kv);
            xmtx = signal_x(kv:kv+buffer-1);
            uv = (reshape(xmtx,[Nh Nv]).')*conj(weights_h(:,kh+1));
            error = signal_d(kv) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1))')*xmtx;
            mu_norm = mu/(norm(uv,2)^2);
            weights_v(:,kv+1) = weights_v(:,kv) + mu_norm*uv*conj(error); 
            weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1));
            mse(kv,1) =  mse(kv,1) + abs(signal_d(kv) - weights'*xmtx).^2;
        end
    end
end
aux3 = toc();
aux3 = aux3/runs;
mse = mse/runs;

txt = ['\mu = 0.5 with K_h = K_v = 500'];
plot(1:Samples, 10*log10(mse),'-','color', [0.9290 0.6940 0.1250], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold off;
title('ATLMS Behavior for Different Learning Steps and Sampling Rate');
xlabel('Iteration');
ylabel('MSE');
xlim([0 10000])
legend_copy = legend("location", "northeast");
set (legend_copy, "fontsize", 12);
grid on;
saveas(gcf,'atlms_sampling.png')

% Time curve

%% TLMS algorithm for benchmark
mu = 0.5;

tic();
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
    end
end
aux0 = toc();
aux0 = aux0/runs;

figure
txt = ['TLMS'];
stem(1,aux0,'-','color', [0.3010 0.7450 0.9330], "linewidth", 3, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['ATLMS Case 1'];
stem(2,aux1,'-','color', [0.4660 0.6740 0.1880], "linewidth", 3, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['ATLMS Case 2'];
stem(3,aux2,'-','color', [0.9290 0.6940 0.1250], "linewidth", 3, "markersize", 8, "DisplayName", txt);
hold on;
txt = ['ATLMS Case 3'];
stem(4,aux3,'-','color', [0.8500 0.3250 0.0980], "linewidth", 3, "markersize", 8, "DisplayName", txt);
hold off;
title('ATLMS Run Time for Different Sampling Rate');
set(gca, 'XTick',[1 2 3 4], 'XTickLabel',{'TLMS','K_h = K_v = 50', 'K_h = K_v = 250', 'K_h = K_v = 500'})
ylabel('Execution Time (s)');
legend_copy = legend("location", "northeast");
set (legend_copy, "fontsize", 10);
grid on;
saveas(gcf,'atlms_time.png')