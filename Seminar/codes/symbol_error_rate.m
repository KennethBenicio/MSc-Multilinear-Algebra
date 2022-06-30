%% Symbol Error Rate (SER) vs. SNR for different filter algorithms
clc;
close all;
clear all;

% System parameters
R = 4;
Nh = 4;
Nv = 4;
Kh = 10;
Kv = 10;
QAM = 4;

mu = 0.1;
runs = 100;
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
    
    % Noise energy
    snr = SNR(ss);
    var_noise = 1/(10^(snr/10));  
    
    for rr = 1:runs
        rr
        % Filter weights
        weights_v = complex(ones(Nv,Samples),ones(Nv,Samples));
        weights_h = complex(ones(Nh,Samples),ones(Nh,Samples));

        % Geometric Channel
        % Angles at the IRS
        theta_AOD_y = 2*pi*randn(1,R)- pi; %AoD y.
        theta_AOD_x = 2*pi*randn(1,R)- pi; %AoD x.
        % Steering vectors AoD x and y
        a_x = 1/sqrt(2) * exp(1i*pi*((0:Nh-1).'*(cos(theta_AOD_x).*sin(theta_AOD_y))));
        a_y = 1/sqrt(2) * exp(1i*pi*((0:Nv-1).'*cos(theta_AOD_y)));
        A   = tensor.mtx_prod_kr(a_y,a_x); % NhNv x R.

        % Received Signal
        signal_d = randi([0,1],Samples,1);
        signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
        signal_x = zeros(Nh*Nv,Samples-R-1);
        for ii = 1:(Samples - R - 1)
            aux = A*signal_d(ii:ii+R-1);
            noise = (sqrt(var_noise/2))*(randn(Nv*Nh,1) + randn(Nv*Nh,1));
            signal_x(:,ii) = aux + noise;
        end

        kh = 1;
        kv = 1;
        % ATLMS algorithm
        for kk = 1:(Kh + Kv):(Kb*(Kh + Kv)) - 5*(Kh + Kv)
            for kh = kk:(kk+Kh-1)
                weights_v(:,kh+1) = weights_v(:,kh);
                xmtx = signal_x(:,kh);
                uh = reshape(xmtx,[Nh Nv])*conj(weights_v(:,kv));
                error = signal_d(kh) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh))')*xmtx;
                mu_norm = mu/(norm(uh,2)^2);
                weights_h(:,kh+1) = weights_h(:,kh) + mu_norm*uh*conj(error); 
                weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh));
                if kh>100
                    ser(ss,1) =  ser(ss,1) + (qamdemod(signal_d(kh),QAM)~=qamdemod(weights'*xmtx,QAM));
                end
            end
            for kv = (kk+Kh):(kk+Kh+Kv-1)
                weights_h(:,kv+1) = weights_h(:,kv);
                xmtx = signal_x(:,kv);
                uv = (reshape(xmtx,[Nh Nv]).')*conj(weights_h(:,kh+1));
                error = signal_d(kv) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1))')*xmtx;
                mu_norm = mu/(norm(uv,2)^2);
                weights_v(:,kv+1) = weights_v(:,kv) + mu_norm*uv*conj(error); 
                weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1));
                if kv>100
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
    
    % Noise energy
    snr = SNR(ss);
    var_noise = 1/(10^(snr/10)); 
    
    for rr = 1:runs
        rr
        % Filter weights
        weights_v = complex(randn(Nv,Samples),randn(Nv,Samples));
        weights_h = complex(randn(Nh,Samples),randn(Nh,Samples));

        % Geometric Channel
        % Angles at the IRS
        theta_AOD_y = 2*pi*randn(1,R)- pi; %AoD y.
        theta_AOD_x = 2*pi*randn(1,R)- pi; %AoD x.
        % Steering vectors AoD x and y
        a_x = 1/sqrt(2) * exp(1i*pi*((0:Nh-1).'*(cos(theta_AOD_x).*sin(theta_AOD_y))));
        a_y = 1/sqrt(2) * exp(1i*pi*((0:Nv-1).'*cos(theta_AOD_y)));
        A   = tensor.mtx_prod_kr(a_y,a_x); % NhNv x R.

        % Received Signal
        signal_d = randi([0,1],Samples,1);
        signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
        signal_x = zeros(Nh*Nv,Samples-R-1);
        for ii = 1:(Samples - R - 1)
            aux = A*signal_d(ii:ii+R-1);
            noise = (sqrt(var_noise/2))*(randn(Nv*Nh,1) + randn(Nv*Nh,1));
            signal_x(:,ii) = aux + noise;
        end

        % TLMS algorithm
        for kk = 1:(Samples - buffer - 1)
            xmtx = signal_x(:,kk);
            uh = reshape(xmtx,[Nh Nv])*conj(weights_v(:,kk));
            uv = (reshape(xmtx,[Nh Nv]).')*conj(weights_h(:,kk));
            error = signal_d(kk) - (tensor.mtx_prod_kron(weights_v(:,kk),weights_h(:,kk))')*xmtx;
            mu_norm = mu/(norm(uh,2)^2 + norm(uv,2)^2);
            weights_h(:,kk+1) = weights_h(:,kk) + mu_norm*uh*conj(error); 
            weights_v(:,kk+1) = weights_v(:,kk) + mu_norm*uv*conj(error); 
            weights = tensor.mtx_prod_kron(weights_v(:,kk),weights_h(:,kk));
            if kk>100
                ser(ss,1) =  ser(ss,1) + (qamdemod(signal_d(kk),QAM)~=qamdemod(weights'*xmtx,QAM));
            end
        end
    end
end
ser = ser/(runs*Samples);

txt = ['TLMS algorithm'];
semilogy(SNR, 10*log10(ser),'-','color', [0.8500 0.3250 0.0980], "linewidth", 3, "markersize", 8, "DisplayName", txt);
hold off;
title('ATLMS and TLMS SER vs. SNR');
xlabel('SNR (dB)');
ylabel('SER (dB)');
legend_copy = legend("location", "northeast");
set (legend_copy, "fontsize", 12);
grid on;
saveas(gcf,'symbol_error_rate.png')
