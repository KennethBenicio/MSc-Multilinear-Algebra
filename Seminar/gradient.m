%% Deterministic Gradient Algorithm
clc;
close all;
clear all;

% Learning rate
mi = 1e-3;
% Filter order
order = 2;
% Number of samples
Samples = 5000;
% Defining the mse error and filter coeficients vectors.
error = zeros(Samples,1);
weights = zeros(order, Samples);

%Wiener Solution
wiener = [0.35;-0.15];

% Defining the energy of the noise vector.
SNR_dB = inf;
SNR_li = 10^(SNR_dB/10);
variance_noise = 1/SNR_li;
noise = sqrt(variance_noise/2).*randn(Samples,1);

% Generating the original signal.
signal_d = randn(Samples,1);

% Convolving the channel and the signal.
Hz = [1 1.6];
signal_x = filter(Hz,1,signal_d);

% Generating the noisy received signal.
signal_x = signal_x + noise;

% Defining the autocorrelation matrix and the cross-correlation vector.
Rx = [3.56, 1.60; 1.60, 3.56;];
p = [1; 0;];

% Obtaining the optimal wiener solution.
wopt = inv(Rx)*p;

% To prevent the missmatch between the filtered signal and the desired
% signal. After some hours of debug I found out that the filtered signal
% was a shifted version of the desired signal.
signal_d = signal_d(order:end,1); 
for ss = 1:(Samples - order - 1)
    % Error between the desired signal and the filtered signal.
    error(ss,1) = signal_d(ss) - weights(:,ss)'*signal_x(ss:ss+order-1);
    % Recursive expression.
    weights(:,ss+1) = weights(:,ss) - 2*mi*(Rx*weights(:,ss) - p);
end

% MSE Curve
figure
semilogy(1:Samples, error.^2,'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 8);
title('Deterministic Gradient Algorithm Behavior');
xlabel('Iterations');
ylabel('MSE');
grid on;
saveas(gcf,'gradient_mse.png')

% Contour
figure
[W0, W1] = meshgrid (-1:0.01:1,-1:0.01:1);
w0 = reshape(W0,[],1);
w1 = reshape(W1,[],1);
[aux,~] = size(w0);
for i = 1:aux
   w = [w0(i); w1(i)]; 
   % We are considering that the desired signal has unitary variance.
   Z(i) = 1 - 2*w.'*p + w.'*Rx*w;
end
Z = reshape(Z,size(W0));
contour(W0,W1,Z);
colormap('gray')
hold on;
for ss = 1:(Samples - order)
   plot(weights(1,ss),weights(2,ss),".-",'color', [0.3010 0.7450 0.9330],"markersize", 8); 
end
plot(wiener(1,1),wiener(2,1),"x",'color', [0.4660 0.6740 0.1880],"markersize", 6); 
hold off;
title('Deterministic Gradient Contour');
xlabel('W_1');
ylabel('W_0');
grid on;
saveas(gcf,'gradient_contour.png')