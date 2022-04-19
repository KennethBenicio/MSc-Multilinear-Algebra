
%% ----- Homework 0 ----- %%
% clc;
% clear;
% close all;
% 
% N = [1 2 3 4 5 6 7 8];
% time1 = zeros(length(N),1);
% time2 = zeros(length(N),1);
% for nn = 1:length(N)
%     for mc = 1:1000
%         A = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
%         B = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
%         tic;
%         inv(tensor.mtx_prod_kron(A,B));
%         aux = toc;
%         time1(nn,1) = time1(nn,1) + aux;
%         tic;
%         tensor.mtx_prod_kron(inv(A),inv(B));
%         aux = toc;
%         time2(nn,1) = time2(nn,1) + aux;
%     end
% end
% time1 = time1/1000;
% time2 = time2/1000;
% 
% figure
% txt = ['(\bf A_{N \times N} \otimes B_{N \times N})^{-1}'];
% semilogy(N,time1,'-d','color', [0 1 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold on;
% txt = ['\bf A^{-1}_{N \times N} \otimes B^{-1}_{N \times N}'];
% semilogy(N,time2,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold off;
% title(['Performance of inverse operation'])
% xlabel('Number of columns')
% ylabel('Time (s)')
% legend_copy = legend("location", "northwest");
% set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
% grid on;
% saveas(gcf,'hw0a1.png')
% 
% N = 2;
% K = [1 2 3 4 5 6 7 8 9 10];
% time1 = zeros(length(K),1);
% time2 = zeros(length(K),1);
% for kk = 1:length(K)
%     for mc = 1:500
%         tic;
%         for ii = 1:K(kk)    
%             if ii == 1
%                 A1 = randn(N,N) + 1j*randn(N,N); 
%                 continue
%             else
%                 A2 = randn(N,N) + 1j*randn(N,N);
%                 A1 = tensor.mtx_prod_kron(A1,A2);
%             end    
%         end
%         inv(A1);
%         aux = toc;
%         time1(kk,1) = time1(kk,1) + aux;
%         
%         tic;
%         for ii = 1:K(kk)    
%             if ii == 1
%                 A1 = randn(N,N) + 1j*randn(N,N); 
%                 continue
%             else
%                 A2 = randn(N,N) + 1j*randn(N,N);
%                 A1 = tensor.mtx_prod_kron(inv(A1),inv(A2));
%             end    
%         end
%         aux = toc;
%         time2(kk,1) = time2(kk,1) + aux;
%     end
% end
% time1 = time1/500;
% time2 = time2/500;
% 
% figure
% txt = ['\bf (\otimes^{K}_{i = 1}A^{(i)}_{4 \times 4})^{-1}'];
% semilogy(K,time1,'-d','color', [0 1 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold on;
% txt = ['\bf \otimes^{K}_{i = 1}(A^{(i)}_{4 \times 4})^{-1}'];
% semilogy(K,time2,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold off;
% title(['Multiple Kronecker products performance'])
% xlabel('Number of products')
% ylabel('Time (s)')
% legend_copy = legend("location", "northwest");
% set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
% grid on;
% saveas(gcf,'hw0a2.png')

%% ----- Homework 1 ----- %%
% clc;
% clear;
% close all;
% 
% N = [2 4 8 16 32 64 128];
% time1 = zeros(length(N),1);
% time2 = zeros(length(N),1);
% for nn = 1:length(N)
%     for mc = 1:1000
%         A = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
%         B = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
%         tic;
%         A.*B;
%         aux = toc;
%         time1(nn,1) = time1(nn,1) + aux;
%         tic;
%         tensor.mtx_prod_had(A,B);
%         aux = toc;
%         time2(nn,1) = time2(nn,1) + aux;
%     end
% end
% time1 = time1/1000;
% time2 = time2/1000;
% 
% figure
% txt = ['Native Function'];
% semilogy(N,time1,'-d','color', [0 1 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold on;
% txt = ['Author Function'];
% semilogy(N,time2,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold off;
% title(['Hadamard product performance'])
% xlabel('Number of columns')
% ylabel('Time (s)')
% legend_copy = legend("location", "southeast");
% set(legend_copy,'Interpreter','tex','location','southeast',"fontsize", 12)
% grid on;
% saveas(gcf,'hw1a1.png')
% 
% N = [2 4 8 16 32 64];
% time1 = zeros(length(N),1);
% time2 = zeros(length(N),1);
% for nn = 1:length(N)
%     for mc = 1:250
%         A = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
%         B = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
%         tic;
%         kron(A,B);
%         aux = toc;
%         time1(nn,1) = time1(nn,1) + aux;
%         tic;
%         tensor.mtx_prod_kron(A,B);
%         aux = toc;
%         time2(nn,1) = time2(nn,1) + aux;
%     end
% end
% time1 = time1/250;
% time2 = time2/250;
% 
% figure
% txt = ['Native Function'];
% semilogy(N,time1,'-d','color', [0 1 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold on;
% txt = ['Author Function'];
% semilogy(N,time2,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold off;
% title(['Kronecker product performance'])
% xlabel('Number of columns')
% ylabel('Time (s)')
% legend_copy = legend("location", "northwest");
% set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
% grid on;
% saveas(gcf,'hw1a2.png')
% 
% N = [2 4 8 16 32 64];
% time1 = zeros(length(N),1);
% time2 = zeros(length(N),1);
% for nn = 1:length(N)
%     for mc = 1:250
%         A = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
%         B = randn(N(nn),N(nn)) + 1j*randn(N(nn),N(nn));
%         tic;
%         kr(A,B);
%         aux = toc;
%         time1(nn,1) = time1(nn,1) + aux;
%         tic;
%         tensor.mtx_prod_kr(A,B);
%         aux = toc;
%         time2(nn,1) = time2(nn,1) + aux;
%     end
% end
% time1 = time1/250;
% time2 = time2/250;
% 
% figure
% txt = ['Native Function'];
% semilogy(N,time1,'-d','color', [0 1 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold on;
% txt = ['Author Function'];
% semilogy(N,time2,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold off;
% title(['Khatri-Rao product performance'])
% xlabel('Number of columns')
% ylabel('Time (s)')
% legend_copy = legend("location", "northwest");
% set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
% grid on;
% saveas(gcf,'hw1a3.png')

%% ----- Homework 2 ----- %%
% clc;
% clear;
% close all;
% 
% R = 2;
% I = [2 4 8 16 32 64];
% time1 = zeros(length(I),1);
% time2 = zeros(length(I),1);
% time3 = zeros(length(I),1);
% for ii = 1:length(I)
%     for mc = 1:250
%         A = randn(I(ii),I(ii)) + 1j*randn(I(ii),I(ii));
%         B = randn(I(ii),I(ii)) + 1j*randn(I(ii),I(ii));
%         tic;
%         pinv(tensor.mtx_prod_kr(A,B));
%         aux = toc;
%         time1(ii,1) = time1(ii,1) + aux;
%         tic;
%         (tensor.mtx_prod_kr(A,B).'*tensor.mtx_prod_kr(A,B))\(tensor.mtx_prod_kr(A,B).');
%         aux = toc;
%         time2(ii,1) = time2(ii,1) + aux;
%         tic;
%         tensor.mtx_prod_had((A.'*A),(B.'*B))\(tensor.mtx_prod_kr(A,B).');
%         aux = toc;
%         time3(ii,1) = time3(ii,1) + aux;
%     end
% end
% time1 = time1/250;
% time2 = time2/250;
% time3 = time3/250;
% 
% figure
% txt = ['1st Method'];
% plot(I,time1,'-d','color', [1 0 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold on;
% txt = ['2nd Method'];
% plot(I,time2,'-d','color', [0 1 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold on;
% txt = ['3rd Method'];
% plot(I,time3,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold off;
% title(['Khatri-Rao product performance considering 2 columns'])
% xlabel('Number of rows')
% ylabel('Time (s)')
% legend_copy = legend("location", "northwest");
% set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
% grid on;
% saveas(gcf,'hw2a1.png')
% 
% R = 4;
% I = [2 4 8 16 32 64];
% time1 = zeros(length(I),1);
% time2 = zeros(length(I),1);
% time3 = zeros(length(I),1);
% for ii = 1:length(I)
%     for mc = 1:250
%         A = randn(I(ii),I(ii)) + 1j*randn(I(ii),I(ii));
%         B = randn(I(ii),I(ii)) + 1j*randn(I(ii),I(ii));
%         tic;
%         pinv(tensor.mtx_prod_kr(A,B));
%         aux = toc;
%         time1(ii,1) = time1(ii,1) + aux;
%         tic;
%         (tensor.mtx_prod_kr(A,B).'*tensor.mtx_prod_kr(A,B))\(tensor.mtx_prod_kr(A,B).');
%         aux = toc;
%         time2(ii,1) = time2(ii,1) + aux;
%         tic;
%         tensor.mtx_prod_had((A.'*A),(B.'*B))\(tensor.mtx_prod_kr(A,B).');
%         aux = toc;
%         time3(ii,1) = time3(ii,1) + aux;
%     end
% end
% time1 = time1/250;
% time2 = time2/250;
% time3 = time3/250;
% 
% figure
% txt = ['1st Method'];
% plot(I,time1,'-d','color', [1 0 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold on;
% txt = ['2nd Method'];
% plot(I,time2,'-d','color', [0 1 0], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold on;
% txt = ['3rd Method'];
% plot(I,time3,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% hold off;
% title(['Khatri-Rao product performance considering 4 columns'])
% xlabel('Number of rows')
% ylabel('Time (s)')
% legend_copy = legend("location", "northwest");
% set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
% grid on;
% saveas(gcf,'hw2a2.png')
% 
% I = 4;
% R = 2;
% N = [2 4 6 8];
% time1 = zeros(length(N),1);
% time2 = zeros(length(N),1);
% for nn = 1:length(N)
%     for mc = 1:250 
%         tic;
%         for ii = 1:N(nn)
%             if ii == 1
%                 A1 = randn(I,R) + 1j*randn(I,R); 
%                 continue
%             else
%                 A2 = randn(I,R) + 1j*randn(I,R);
%                 A1 = tensor.mtx_prod_kron(A1,A2);
%             end    
%         end
%         aux = toc;
%         time1(nn,1) = time1(nn,1) + aux;
%     end
% end
% time1 = time1/250;
% 
% figure
% txt = ['\bf ⋄^{N}_{n = 1} (A^{(n)}_{4 \times 2})'];
% semilogy(N,time1,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% title(['Multiple Khatri-Rao products performance'])
% xlabel('Number of Products')
% ylabel('Time (s)')
% legend_copy = legend("location", "northwest");
% set(legend_copy,'Interpreter','tex','location','northwest',"fontsize", 12)
% grid on;
% saveas(gcf,'hw2a3.png')

%% ----- Homework 3 ----- %%
% clc;
% clear;
% close all;
% 
% A = randn(4,2) + 1j*randn(4,2);
% B = randn(6,2) + 1j*randn(6,2);
% X = tensor.mtx_prod_kr(A,B);
% [Ahat,Bhat] = tensor.LSKRF(X,4,6);
% Xhat = tensor.mtx_prod_kr(Ahat,Bhat);
% 
% nmsex = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
% nmsex = 20*log10(nmsex)
% nmsea = (norm(A- Ahat,'fro')^2)/(norm(A,'fro')^2);
% nmsea = 20*log10(nmsea)
% nmseb = (norm(B- Bhat,'fro')^2)/(norm(B,'fro')^2);
% nmseb = 20*log10(nmseb)
% 
% I = 10;
% J = 10;
% R =  4;
% SNR = [0 5 10 15 20 25 30];
% nmse = zeros(length(SNR),1);
% for snr = 1:length(SNR)
%     for mc = 1:1000
%         var_noise = 1/(10^(SNR(snr)/10));
%         noise = sqrt(var_noise/2)*(randn(I*J,R) + 1j*randn(I*J,R));
%         
%         A = randn(I,R) + 1j*randn(I,R);
%         B = randn(J,R) + 1j*randn(J,R);
%         X = tensor.mtx_prod_kr(A,B);
%         X = X + noise;
%         
%         [Ahat,Bhat] = tensor.LSKRF(X,I,J);
%         Xhat = tensor.mtx_prod_kr(Ahat,Bhat);
%         aux = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
%         nmse(snr,1) = nmse(snr,1) + 20*log10(aux);
%     end
% end
% nmse  = nmse/1000;
% 
% figure
% txt = ['I = ' num2str(I), ', J = ' num2str(J), ' and R = ' num2str(R)];
% plot(SNR,nmse,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% title(['LSKRF performance under imperfect scenario'])
% xlabel('SNR (dB)')
% ylabel('NMSE (dB)')
% legend_copy = legend("location", "northwest");
% set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
% grid on;
% saveas(gcf,'hw3a1.png')
% 
% I = 30;
% J = 10;
% R =  4;
% SNR = [0 5 10 15 20 25 30];
% nmse = zeros(length(SNR),1);
% for snr = 1:length(SNR)
%     for mc = 1:1000
%         var_noise = 1/(10^(SNR(snr)/10));
%         noise = sqrt(var_noise/2)*(randn(I*J,R) + 1j*randn(I*J,R));
%         
%         A = randn(I,R) + 1j*randn(I,R);
%         B = randn(J,R) + 1j*randn(J,R);
%         X = tensor.mtx_prod_kr(A,B);
%         X = X + noise;
%         
%         [Ahat,Bhat] = tensor.LSKRF(X,I,J);
%         Xhat = tensor.mtx_prod_kr(Ahat,Bhat);
%         aux = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
%         nmse(snr,1) = nmse(snr,1) + 20*log10(aux);
%     end
% end
% nmse  = nmse/1000;
% 
% figure
% txt = ['I = ' num2str(I), ', J = ' num2str(J), ' and R = ' num2str(R)];
% plot(SNR,nmse,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% title(['LSKRF performance under imperfect scenario'])
% xlabel('SNR (dB)')
% ylabel('NMSE (dB)')
% legend_copy = legend("location", "northwest");
% set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
% grid on;
% saveas(gcf,'hw3a2.png')

%% ----- Homework 4 ----- %%
% clc;
% clear;
% close all;
% 
% A = randn(4,2) + 1j*randn(4,2);
% B = randn(6,3) + 1j*randn(6,3);
% X = tensor.mtx_prod_kron(A,B);
% [Ahat,Bhat] = tensor.LSKronF(X,4,2,6,3);
% Xhat = tensor.mtx_prod_kron(Ahat,Bhat);
% 
% nmsex = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
% nmsex = 20*log10(nmsex)
% nmsea = (norm(A- Ahat,'fro')^2)/(norm(A,'fro')^2);
% nmsea = 20*log10(nmsea)
% nmseb = (norm(B- Bhat,'fro')^2)/(norm(B,'fro')^2);
% nmseb = 20*log10(nmseb)
% 
% I = 2;
% J = 4;
% P = 3;
% Q = 5;
% SNR = [0 5 10 15 20 25 30];
% nmse = zeros(length(SNR),1);
% for snr = 1:length(SNR)
%     for mc = 1:1000
%         var_noise = 1/(10^(SNR(snr)/10));
%         noise = sqrt(var_noise/2)*(randn(I*J,P*Q) + 1j*randn(I*J,P*Q));
%         
%         A = randn(I,P) + 1j*randn(I,P);
%         B = randn(J,Q) + 1j*randn(J,Q);
%         X = tensor.mtx_prod_kron(A,B);
%         X = X + noise;
%         
%         [Ahat,Bhat] = tensor.LSKronF(X,I,P,J,Q);
%         Xhat = tensor.mtx_prod_kron(Ahat,Bhat);
%         aux = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
%         nmse(snr,1) = nmse(snr,1) + 20*log10(aux);
%     end
% end
% nmse  = nmse/1000;
% 
% figure
% txt = ['I = ' num2str(I), ', J = ' num2str(J), ', P = ' num2str(P), ' and Q = ' num2str(Q)];
% plot(SNR,nmse,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% title(['LSKronF performance under imperfect scenario'])
% xlabel('SNR (dB)')
% ylabel('NMSE (dB)')
% legend_copy = legend("location", "northwest");
% set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
% grid on;
% saveas(gcf,'hw4a1.png')
% 
% I = 4;
% J = 8;
% P = 3;
% Q = 5;
% SNR = [0 5 10 15 20 25 30];
% nmse = zeros(length(SNR),1);
% for snr = 1:length(SNR)
%     for mc = 1:1000
%         var_noise = 1/(10^(SNR(snr)/10));
%         noise = sqrt(var_noise/2)*(randn(I*J,P*Q) + 1j*randn(I*J,P*Q));
%         
%         A = randn(I,P) + 1j*randn(I,P);
%         B = randn(J,Q) + 1j*randn(J,Q);
%         X = tensor.mtx_prod_kron(A,B);
%         X = X + noise;
%         
%         [Ahat,Bhat] = tensor.LSKronF(X,I,P,J,Q);
%         Xhat = tensor.mtx_prod_kron(Ahat,Bhat);
%         aux = (norm(X- Xhat,'fro')^2)/(norm(X,'fro')^2);
%         nmse(snr,1) = nmse(snr,1) + 20*log10(aux);
%     end
% end
% nmse  = nmse/1000;
% 
% figure
% txt = ['I = ' num2str(I), ', J = ' num2str(J), ', P = ' num2str(P), ' and Q = ' num2str(Q)];
% plot(SNR,nmse,'-d','color', [0 0 1], "linewidth", 2, "markersize", 8, "DisplayName", txt);
% title(['LSKronF performance under imperfect scenario'])
% xlabel('SNR (dB)')
% ylabel('NMSE (dB)')
% legend_copy = legend("location", "northwest");
% set(legend_copy,'Interpreter','tex','location','northeast',"fontsize", 12)
% grid on;
% saveas(gcf,'hw4a2.png')

%% ----- Homework 5 ----- %%
clc;
clear;
close all;

% My tensor for testing
disp('Generic Tensor:');
tsr_xxx = load('test_unfold_fold.mat').xxx
dimension = size(tsr_xxx);

% Unfolding
disp('1st mode unfolding:');
tsr_xxx_1 = tensor.unfold(tsr_xxx,1)
disp('2nd mode unfolding:');
tsr_xxx_2 = tensor.unfold(tsr_xxx,2)
disp('3rd mode unfolding:');
tsr_xxx_3 = tensor.unfold(tsr_xxx,3)

% Folding
disp('1st mode folding:');
tsr_xxx_1 = tensor.fold(tsr_xxx_1,dimension,1)
disp('2nd mode folding:');
tsr_xxx_2 = tensor.fold(tsr_xxx_2,dimension,2)
disp('3rd mode folding:');
tsr_xxx_3 = tensor.fold(tsr_xxx_3,dimension,3)

%% ----- Homework 6 ----- %%

%% ----- Homework 7 ----- %%

%% ----- Homework 8 ----- %%

%% ----- Homework 9 ----- %%