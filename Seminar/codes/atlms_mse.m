%% Alternating Tensor LMS Algorithm
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
runs = 100;
buffer = Nh*Nv;
Samples = 15000;

% Noise energy
snr = 30;
var_noise = 1/(10^(snr/10));  

%% 0.1 Learning Step
mu = 0.1;

% Monter Carlo simulations
mse = zeros(Samples,1);
Kb = floor(Samples/(Kh + Kv));
for rr = 1:runs
	rr
    % Filter weights
    weights_v = complex(ones(Nv,Samples),ones(Nv,Samples));
    weights_h = complex(ones(Nh,Samples),ones(Nh,Samples));
    
    % Geometric Channel
    % Angles at the IRS
    theta_AOD_y = pi*rand(1,R)- pi/2; %AoD y.
    theta_AOD_x = pi*rand(1,R)- pi/2; %AoD x.
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
            mse(kh,1) =  mse(kh,1) + abs(signal_d(kh) - weights'*xmtx).^2;
        end
        for kv = (kk+Kh):(kk+Kh+Kv-1)
            weights_h(:,kv+1) = weights_h(:,kv);
            xmtx = signal_x(:,kv);
            uv = (reshape(xmtx,[Nh Nv]).')*conj(weights_h(:,kh+1));
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
txt = ['\mu = 0.1'];
plot(1:Samples, 10*log10(mse),'-','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;

% %% 0.5 Learning Step
% mu = 0.5;
% 
% % Monter Carlo simulations
% mse = zeros(Samples,1);
% Kb = floor(Samples/(Kh + Kv));
% for rr = 1:runs
% 	rr
%     % Filter weights
%     weights_v = complex(ones(Nv,Samples),ones(Nv,Samples));
%     weights_h = complex(ones(Nh,Samples),ones(Nh,Samples));
%     
%     % Geometric Channel
%     % Angles at the IRS
%     theta_AOD_y = pi*rand(1,R)- pi/2; %AoD y.
%     theta_AOD_x = pi*rand(1,R)- pi/2; %AoD x.
%     % Steering vectors AoD x and y
%     a_x = 1/sqrt(2) * exp(1i*pi*((0:Nh-1).'*(cos(theta_AOD_x).*sin(theta_AOD_y))));
%     a_y = 1/sqrt(2) * exp(1i*pi*((0:Nv-1).'*cos(theta_AOD_y)));
%     A   = tensor.mtx_prod_kr(a_y,a_x); % NhNv x R.
%     
%     % Received Signal
%     signal_d = randi([0,1],Samples,1);
%     signal_d = 1/sqrt(2) * qammod(signal_d,QAM);
%     signal_x = zeros(Nh*Nv,Samples-R-1);
%     for ii = 1:(Samples - R - 1)
%         aux = A*signal_d(ii:ii+R-1);
%         noise = (sqrt(var_noise/2))*(randn(Nv*Nh,1) + randn(Nv*Nh,1));
%         signal_x(:,ii) = aux + noise;
%     end
%     
%     kh = 1;
%     kv = 1;
%     % ATLMS algorithm
%     for kk = 1:(Kh + Kv):(Kb*(Kh + Kv)) - 5*(Kh + Kv)
%         for kh = kk:(kk+Kh-1)
%             weights_v(:,kh+1) = weights_v(:,kh);
%             xmtx = signal_x(:,kh);
%             uh = reshape(xmtx,[Nh Nv])*conj(weights_v(:,kv));
%             error = signal_d(kh) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh))')*xmtx;
%             mu_norm = mu/(norm(uh,2)^2);
%             weights_h(:,kh+1) = weights_h(:,kh) + mu_norm*uh*conj(error); 
%             weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh));
%             mse(kh,1) =  mse(kh,1) + abs(signal_d(kh) - weights'*xmtx).^2;
%         end
%         for kv = (kk+Kh):(kk+Kh+Kv-1)
%             weights_h(:,kv+1) = weights_h(:,kv);
%             xmtx = signal_x(:,kv);
%             uv = (reshape(xmtx,[Nh Nv]).')*conj(weights_h(:,kh+1));
%             error = signal_d(kv) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1))')*xmtx;
%             mu_norm = mu/(norm(uv,2)^2);
%             weights_v(:,kv+1) = weights_v(:,kv) + mu_norm*uv*conj(error); 
%             weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1));
%             mse(kv,1) =  mse(kv,1) + abs(signal_d(kv) - weights'*xmtx).^2;
%         end
%     end
% end
% mse = mse/runs;
% 
% txt = ['\mu = 0.5'];
% plot(1:Samples, 10*log10(mse),'-','color', [0.8500 0.3250 0.0980], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold on;

%% 0.05 Learning Step
mu = 0.05;

% Monter Carlo simulations
mse = zeros(Samples,1);
Kb = floor(Samples/(Kh + Kv));
for rr = 1:runs
	rr
    % Filter weights
    weights_v = complex(ones(Nv,Samples),ones(Nv,Samples));
    weights_h = complex(ones(Nh,Samples),ones(Nh,Samples));
    
    % Geometric Channel
    % Angles at the IRS
    theta_AOD_y = pi*rand(1,R)- pi/2; %AoD y.
    theta_AOD_x = pi*rand(1,R)- pi/2; %AoD x.
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
            mse(kh,1) =  mse(kh,1) + abs(signal_d(kh) - weights'*xmtx).^2;
        end
        for kv = (kk+Kh):(kk+Kh+Kv-1)
            weights_h(:,kv+1) = weights_h(:,kv);
            xmtx = signal_x(:,kv);
            uv = (reshape(xmtx,[Nh Nv]).')*conj(weights_h(:,kh+1));
            error = signal_d(kv) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1))')*xmtx;
            mu_norm = mu/(norm(uv,2)^2);
            weights_v(:,kv+1) = weights_v(:,kv) + mu_norm*uv*conj(error); 
            weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1));
            mse(kv,1) =  mse(kv,1) + abs(signal_d(kv) - weights'*xmtx).^2;
        end
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
Kb = floor(Samples/(Kh + Kv));
for rr = 1:runs
	rr
    % Filter weights
    weights_v = complex(ones(Nv,Samples),ones(Nv,Samples));
    weights_h = complex(ones(Nh,Samples),ones(Nh,Samples));
    
    % Geometric Channel
    % Angles at the IRS
    theta_AOD_y = pi*rand(1,R)- pi/2; %AoD y.
    theta_AOD_x = pi*rand(1,R)- pi/2; %AoD x.
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
            mse(kh,1) =  mse(kh,1) + abs(signal_d(kh) - weights'*xmtx).^2;
        end
        for kv = (kk+Kh):(kk+Kh+Kv-1)
            weights_h(:,kv+1) = weights_h(:,kv);
            xmtx = signal_x(:,kv);
            uv = (reshape(xmtx,[Nh Nv]).')*conj(weights_h(:,kh+1));
            error = signal_d(kv) - (tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1))')*xmtx;
            mu_norm = mu/(norm(uv,2)^2);
            weights_v(:,kv+1) = weights_v(:,kv) + mu_norm*uv*conj(error); 
            weights = tensor.mtx_prod_kron(weights_v(:,kv),weights_h(:,kh+1));
            mse(kv,1) =  mse(kv,1) + abs(signal_d(kv) - weights'*xmtx).^2;
        end
    end
end
mse = mse/runs;

txt = ['\mu = 0.005'];
plot(1:Samples, 10*log10(mse),'-','color', [0.9290 0.6940 0.1250], "linewidth", 2, "markersize", 8, "DisplayName", txt);
hold on;
title('ATLMS Behavior for Different Learning Steps');
xlabel('Iteration');
ylabel('MSE');
legend_copy = legend("location", "northeast");
set (legend_copy, "fontsize", 12);
grid on;
saveas(gcf,'atlms_mse.png')