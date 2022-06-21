%% Symbol Error Rate (SER) vs. SNR for different filter algorithms
clc;
close all;
clear all;

% System parameters
Nh = 5;
Nv = 5;
Kh = 5;
Kv = 5;
order = 10;

QAM = 2;
mu = 0.1;
runs = 5000;
buffer = Nh*Nv;
Samples = 10000;

% SNR vector
SNR = [0 5 10 15 20 25 30];

%% ATLMS SER

% Monter Carlo simulations
ser = zeros(length(SNR),1);
Kb = floor(Samples/(Kh + Kv));
for ss = 1:length(SNR)
    ss
    for rr = 1:runs
        rr
        % Filter weights
        weights_v = complex(ones(Nv,Samples),ones(Nv,Samples));
        weights_h = complex(ones(Nh,Samples),ones(Nh,Samples));

        % M-QAM source signal
        signal_d = randi([0,1],Samples,1);
        signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
        % Adding noise to the M-QAM source signal
        var_noise = 1/(10^(SNR(ss)/10));  
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
                if kh>250
                    ser(ss,1) =  ser(ss,1) + (qamdemod(signal_d(kh),QAM)~=qamdemod(weights'*xmtx,QAM));
                end
            end
            for kv = (kk+Kh):(kk+Kh+Kv-1)
                weights_h(:,kv+1) = weights_h(:,kv);
                xmtx = signal_x(kv:kv+buffer-1);
                uv = (reshape(xmtx,[Nh Nv]).')*conj(weights_h(:,kh+1));
                error = signal_d(kv) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1))')*xmtx;
                mu_norm = mu/(norm(uv,2)^2);
                weights_v(:,kv+1) = weights_v(:,kv) + mu_norm*uv*conj(error); 
                weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1));
                if kv>250
                    ser(ss,1) =  ser(ss,1) + (qamdemod(signal_d(kv),QAM)~=qamdemod(weights'*xmtx,QAM));
                end
            end
        end
    end
end
ser = ser/(runs*Samples);

% SER Curve
figure
txt = ['ATLMS algorithm'];
semilogy(SNR, 10*log10(ser),'-','color', [0.3010 0.7450 0.9330], "linewidth", 3, "markersize", 8, "DisplayName", txt);
hold on;

%% TLMS SER

% Monter Carlo simulations
ser = zeros(length(SNR),1);
for ss = 1:length(SNR)
    ss
    for rr = 1:runs
        rr
        % Filter weights
        weights_v = complex(ones(Nv,Samples),ones(Nv,Samples));
        weights_h = complex(ones(Nh,Samples),ones(Nh,Samples));

        % M-QAM source signal
        signal_d = randi([0,1],Samples,1);
        signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
        % Adding noise to the M-QAM source signal
        var_noise = 1/(10^(SNR(ss)/10));  
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
            if kk>250
                ser(ss,1) =  ser(ss,1) + (qamdemod(signal_d(kk),QAM)~=qamdemod(weights'*xmtx,QAM));
            end
        end
    end
end
ser = ser/(runs*Samples);

txt = ['TLMS algorithm'];
semilogy(SNR, 10*log10(ser),'-','color', [0.8500 0.3250 0.0980], "linewidth", 3, "markersize", 8, "DisplayName", txt);
hold on;

%% LMS SER

% Monter Carlo simulations
ser = zeros(length(SNR),1);
for ss = 1:length(SNR)
    ss
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
        
        % LMS algorithm
        for kk = 1:(Samples - order - 1)
            % Error between the desired signal and the filtered signal.
            error = signal_d(kk) - weights(:,kk)' * signal_x(kk:kk+order-1); 
            % Recursive expression.
            weights(:,kk+1) = weights(:,kk) + 2 * mu * error * signal_x(kk:kk+order-1);
            if kk>250
                ser(ss,1) =  ser(ss,1) + (qamdemod(signal_d(kk),QAM)~=qamdemod(weights(:,kk+1)'*signal_x(kk:kk+order-1),QAM));
            end
        end
    end
end
ser = ser/(runs*Samples);

txt = ['LMS algorithm'];
semilogy(SNR, 10*log10(ser),'-','color', [0.4660 0.6740 0.1880], "linewidth", 3, "markersize", 8, "DisplayName", txt);
hold on;

%% NLMS SER

% Monter Carlo simulations
ser = zeros(length(SNR),1);
for ss = 1:length(SNR)
    ss
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
            if kk>250
                ser(ss,1) =  ser(ss,1) + (qamdemod(signal_d(kk),QAM)~=qamdemod(weights(:,kk+1)'*signal_x(kk:kk+order-1),QAM));
            end
        end
    end
end
ser = ser/(runs*Samples);

txt = ['NLMS algorithm'];
semilogy(SNR, 10*log10(ser),'-','color', [0.9290 0.6940 0.1250], "linewidth", 3, "markersize", 8, "DisplayName", txt);
hold off;
title('ATLMS SER vs. SNR');
xlabel('SNR (dB)');
ylabel('SER (dB)');
legend_copy = legend("location", "northeast");
set (legend_copy, "fontsize", 12);
grid on;
saveas(gcf,'symbol_error_rate.png')