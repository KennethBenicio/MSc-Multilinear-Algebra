%% Normalized LMS Algorithm
clc;
close all;
clear all;

% System parameters
QAM = 2;
runs = 2500;
order = 1;
Samples = 15000;

% Noise energy
snr = 30;
var_noise = 1/(10^(snr/10));   

%% 0.1 Learning Step
mu = 0.1;

% Monter Carlo simulations
mse = zeros(Samples,1);
for rr = 1:runs
	rr
    % Filter weights
    weights = complex(zeros(order, Samples),0);
    
    % M-QAM source signal
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
    
    % Adding noise to the M-QAM source signal
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    
    % NLMS algorithm
    for kk = 1:(Samples - order - 1)
        mu_normalized = mu/(norm(signal_x(kk:kk+order-1)));
       % Error between the desired signal and the filtered signal.
        error(kk) = signal_d(kk) - weights(:,kk)' * signal_x(kk:kk+order-1); 
        % Recursive expression.
        weights(:,kk+1) = weights(:,kk) + mu_normalized * error(kk) * signal_x(kk:kk+order-1);
        mse(kk,1) =  mse(kk,1) + abs(signal_d(kk) - weights(:,kk+1)'*signal_x(kk:kk+order-1)).^2;
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

% Monter Carlo simulations
mse = zeros(Samples,1);
for rr = 1:runs
    % Filter weights
    weights = complex(zeros(order, Samples),0);
    
    % M-QAM source signal
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
    
    % Adding noise to the M-QAM source signal
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    
    % NLMS algorithm
    for kk = 1:(Samples - order - 1)
        mu_normalized = mu/(norm(signal_x(kk:kk+order-1)));
       % Error between the desired signal and the filtered signal.
        error(kk) = signal_d(kk) - weights(:,kk)' * signal_x(kk:kk+order-1); 
        % Recursive expression.
        weights(:,kk+1) = weights(:,kk) + mu_normalized * error(kk) * signal_x(kk:kk+order-1);
        mse(kk,1) =  mse(kk,1) + abs(signal_d(kk) - weights(:,kk+1)'*signal_x(kk:kk+order-1)).^2;
    end
end
mse = mse/runs;

txt = ['\mu = 0.5'];
plot(1:Samples, 10*log10(mse),'-','color', [0.8500 0.3250 0.0980], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;

%% 0.05 Learning Step
mu = 0.05;

% Monter Carlo simulations
mse = zeros(Samples,1);
for rr = 1:runs
    % Filter weights
    weights = complex(zeros(order, Samples),0);
    
    % M-QAM source signal
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
    
    % Adding noise to the M-QAM source signal
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    
    % NLMS algorithm
    for kk = 1:(Samples - order - 1)
        mu_normalized = mu/(norm(signal_x(kk:kk+order-1)));
       % Error between the desired signal and the filtered signal.
        error(kk) = signal_d(kk) - weights(:,kk)' * signal_x(kk:kk+order-1); 
        % Recursive expression.
        weights(:,kk+1) = weights(:,kk) + mu_normalized * error(kk) * signal_x(kk:kk+order-1);
        mse(kk,1) =  mse(kk,1) + abs(signal_d(kk) - weights(:,kk+1)'*signal_x(kk:kk+order-1)).^2;
    end
end
mse = mse/runs;

txt = ['\mu = 0.05'];
plot(1:Samples, 10*log10(mse),'-','color', [0.4660 0.6740 0.1880], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;

%% 0.005 Learning Step
mu = 0.005;

% Monter Carlo simulations
mse = zeros(Samples,1);
for rr = 1:runs
    % Filter weights
    weights = complex(zeros(order, Samples),0);
    
    % M-QAM source signal
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
    
    % Adding noise to the M-QAM source signal
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
    
    % NLMS algorithm
    for kk = 1:(Samples - order - 1)
        mu_normalized = mu/(norm(signal_x(kk:kk+order-1)));
       % Error between the desired signal and the filtered signal.
        error(kk) = signal_d(kk) - weights(:,kk)' * signal_x(kk:kk+order-1); 
        % Recursive expression.
        weights(:,kk+1) = weights(:,kk) + mu_normalized * error(kk) * signal_x(kk:kk+order-1);
        mse(kk,1) =  mse(kk,1) + abs(signal_d(kk) - weights(:,kk+1)'*signal_x(kk:kk+order-1)).^2;
    end
end
mse = mse/runs;

txt = ['\mu = 0.005'];
plot(1:Samples, 10*log10(mse),'-','color', [0.9290 0.6940 0.1250], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
title('NLMS Behavior for Different Learning Steps');
xlabel('Iteration');
ylabel('MSE');
legend_copy = legend("location", "northeast");
set (legend_copy, "fontsize", 12);
grid on;
saveas(gcf,'nlms_mse.png')
