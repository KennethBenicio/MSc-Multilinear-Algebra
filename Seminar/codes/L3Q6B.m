clc;
close all;
clear all;

% Learning rate
mi = 1e-3;
% Filter order
order = 15;

% If we want to evaluate the impact of the training size then we should
% apply the same transmitted sequence to all the cases

% Number of samples
Samples = 5000;
% Defining the energy of the noise vector.
SNR = 30;
QAM = 16;
signal_d = randi([0,QAM - 1],[Samples 1]); % The same pilot for every pilot frame and block.
signal_d = (1/sqrt(2)) * qammod(signal_d,QAM); % 4-QAM Pilot Signal.

% Convolving the channel and the signal.
Hz = [0.5 1.2 1.5 -1];
signal_x = filter(Hz,1,signal_d);

% Training noise
snr = 10^(SNR/10);
energy_symbol = mean(abs(signal_x(:)).^2); % Energy symbol pilot. 
var_noise = energy_symbol .*  1/snr; % Variance of the noise.
noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

% Generating the noisy received signal.
signal_x = signal_x + noise;

%% Training with 50 Samples
% Number of samples
Samples = 50;
% Defining the mse error and filter coeficients vectors.
error = zeros(Samples,1);
weights = zeros(order, Samples);

% Defining the energy of the noise vector.
SNR = 30;
QAM_train = 4;
signal_d_train = randi([0,QAM_train - 1],[Samples 1]); % The same pilot for every pilot frame and block.
signal_d_train = (1/sqrt(2)) * qammod(signal_d_train,QAM_train); % 4-QAM Pilot Signal.

% Convolving the channel and the signal.
Hz = [0.5 1.2 1.5 -1];
signal_x_train = filter(Hz,1,signal_d_train);

% Training noise
snr = 10^(SNR/10);
energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
var_noise = energy_symbol .*  1/snr; % Variance of the noise.
noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

% Generating the noisy received signal.
signal_x_train = signal_x_train + noise;

for s = order:Samples
    aux = signal_x_train(s:-1:s-order+1);
    error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
    % Recursive expression.
    weights(:,s+1) = weights(:,s) + 2* mi * conj(error(s)) * aux;
end

% Number of samples
Samples = 5000;
% Defining the mse error and filter coeficients vectors.
error = zeros(Samples,1);
aux = weights(:,s);
weights = zeros(order, Samples);
weights(:,1) = aux;

signal_d_hat_50 = zeros(size(signal_d));
for s = order:Samples
    aux = signal_x(s:-1:s-order+1);
    signal_d_hat_50(s) = weights(:,s)'*aux;
    error(s) = signal_d(s-order+1) - signal_d_hat_50(s);
    % Recursive expression.
    weights(:,s+1) = weights(:,s) + 2 * mi * conj(error(s)) * aux;
end

%% Training with 150 Samples
% Number of samples
Samples = 150;
% Defining the mse error and filter coeficients vectors.
error = zeros(Samples,1);
weights = zeros(order, Samples);

% Defining the energy of the noise vector.
SNR = 30;
QAM_train = 4;
signal_d_train = randi([0,QAM_train - 1],[Samples 1]); % The same pilot for every pilot frame and block.
signal_d_train = (1/sqrt(2)) * qammod(signal_d_train,QAM_train); % 4-QAM Pilot Signal.

% Convolving the channel and the signal.
Hz = [0.5 1.2 1.5 -1];
signal_x_train = filter(Hz,1,signal_d_train);

% Training noise
snr = 10^(SNR/10);
energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
var_noise = energy_symbol .*  1/snr; % Variance of the noise.
noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

% Generating the noisy received signal.
signal_x_train = signal_x_train + noise;

for s = order:Samples
    aux = signal_x_train(s:-1:s-order+1);
    error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
    % Recursive expression.
    weights(:,s+1) = weights(:,s) + 2* mi * conj(error(s)) * aux;
end

% Number of samples
Samples = 5000;
% Defining the mse error and filter coeficients vectors.
error = zeros(Samples,1);
aux = weights(:,s);
weights = zeros(order, Samples);
weights(:,1) = aux;

signal_d_hat_150 = zeros(size(signal_d));
for s = order:Samples
    aux = signal_x(s:-1:s-order+1);
    signal_d_hat_150(s) = weights(:,s)'*aux;
    error(s) = signal_d(s-order+1) - signal_d_hat_150(s);
    % Recursive expression.
    weights(:,s+1) = weights(:,s) + 2 * mi * conj(error(s)) * aux;
end

%% Training with 350 Samples
% Number of samples
Samples = 350;
% Defining the mse error and filter coeficients vectors.
error = zeros(Samples,1);
weights = zeros(order, Samples);

% Defining the energy of the noise vector.
SNR = 30;
QAM_train = 4;
signal_d_train = randi([0,QAM_train - 1],[Samples 1]); % The same pilot for every pilot frame and block.
signal_d_train = (1/sqrt(2)) * qammod(signal_d_train,QAM_train); % 4-QAM Pilot Signal.

% Convolving the channel and the signal.
Hz = [0.5 1.2 1.5 -1];
signal_x_train = filter(Hz,1,signal_d_train);

% Training noise
snr = 10^(SNR/10);
energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
var_noise = energy_symbol .*  1/snr; % Variance of the noise.
noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

% Generating the noisy received signal.
signal_x_train = signal_x_train + noise;

for s = order:Samples
    aux = signal_x_train(s:-1:s-order+1);
    error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
    % Recursive expression.
    weights(:,s+1) = weights(:,s) + 2* mi * conj(error(s)) * aux;
end

% Number of samples
Samples = 5000;
% Defining the mse error and filter coeficients vectors.
error = zeros(Samples,1);
aux = weights(:,s);
weights = zeros(order, Samples);
weights(:,1) = aux;

signal_d_hat_300 = zeros(size(signal_d));
for s = order:Samples
    aux = signal_x(s:-1:s-order+1);
    signal_d_hat_300(s) = weights(:,s)'*aux;
    error(s) = signal_d(s-order+1) - signal_d_hat_300(s);
    % Recursive expression.
    weights(:,s+1) = weights(:,s) + 2 * mi * conj(error(s)) * aux;
end

%% Training with 500 Samples
% Number of samples
Samples = 500;
% Defining the mse error and filter coeficients vectors.
error = zeros(Samples,1);
weights = zeros(order, Samples);

% Defining the energy of the noise vector.
SNR = 30;
QAM_train = 4;
signal_d_train = randi([0,QAM_train - 1],[Samples 1]); % The same pilot for every pilot frame and block.
signal_d_train = (1/sqrt(2)) * qammod(signal_d_train,QAM_train); % 4-QAM Pilot Signal.

% Convolving the channel and the signal.
Hz = [0.5 1.2 1.5 -1];
signal_x_train = filter(Hz,1,signal_d_train);

% Training noise
snr = 10^(SNR/10);
energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
var_noise = energy_symbol .*  1/snr; % Variance of the noise.
noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

% Generating the noisy received signal.
signal_x_train = signal_x_train + noise;

for s = order:Samples
    aux = signal_x_train(s:-1:s-order+1);
    error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
    % Recursive expression.
    weights(:,s+1) = weights(:,s) + 2* mi * conj(error(s)) * aux;
end

% Number of samples
Samples = 5000;
% Defining the mse error and filter coeficients vectors.
error = zeros(Samples,1);
aux = weights(:,s);
weights = zeros(order, Samples);
weights(:,1) = aux;

signal_d_hat_500 = zeros(size(signal_d));
for s = order:Samples
    aux = signal_x(s:-1:s-order+1);
    signal_d_hat_500(s) = weights(:,s)'*aux;
    error(s) = signal_d(s-order+1) - signal_d_hat_500(s);
    % Recursive expression.
    weights(:,s+1) = weights(:,s) + 2 * mi * conj(error(s)) * aux;
end

%% Temporal Evolution
% https://www.mathworks.com/help/comm/gs/examine-16-qam-using-matlab.html
figure
subplot(221)
plot(signal_d_hat_50,'.','color', [0.3010 0.7450 0.9330],"markersize", 3);
hold on;
plot(signal_d,"*",'color', [1 0 0],"markersize", 9)
hold off;
xlim([-3 3]);
ylim([-3 3]);
title('Training with 50 Samples');
xlabel('In Phase');
ylabel('Quadrature');
grid on;
subplot(222)
plot(signal_d_hat_150,'.','color', [0.3010 0.7450 0.9330],"markersize", 3);
hold on;
plot(signal_d,"*",'color', [1 0 0],"markersize", 9)
hold off;
xlim([-3 3]);
ylim([-3 3]);
title('Training with 150 Samples');
xlabel('In Phase');
ylabel('Quadrature');
grid on;
subplot(223)
plot(signal_d_hat_300,'.','color', [0.3010 0.7450 0.9330],"markersize", 3);
hold on;
plot(signal_d,"*",'color', [1 0 0],"markersize", 9)
hold off;
xlim([-3 3]);
ylim([-3 3]);
title('Training with 300 Samples');
xlabel('In Phase');
ylabel('Quadrature');
grid on;
subplot(224)
plot(signal_d_hat_500,'.','color', [0.3010 0.7450 0.9330],"markersize", 3);
hold on;
plot(signal_d,"*",'color', [1 0 0],"markersize", 9)
hold off;
xlim([-3 3]);
ylim([-3 3]);
title('Training with 500 Samples');
xlabel('In Phase');
ylabel('Quadrature');
grid on;
saveas(gcf,'L3Q6_B_t.png')