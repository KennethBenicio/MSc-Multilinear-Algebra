%% Normalized LMS Algorithm
clc;
close all;
clear all;

runs = 5;
mu = 0.005;

QAM = 2;
snr = 30;
var_noise = 1/(10^(snr/10));     

order = 10;
Samples = 10000;
mse = zeros(Samples,1);
weights = complex(zeros(order, Samples),0);
for rr = 1:runs
    signal_d = randi([0,1],Samples,1);
    signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
    noise = (sqrt(var_noise/2))*(randn(Samples,1) + randn(Samples,1));
    signal_x = signal_d + noise;
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
plot(1:Samples, 10*log10(mse),'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 8);
title('NLMS Behavior');
xlabel('Samples');
ylabel('MSE');
grid on;
saveas(gcf,'nlms_mse.png')