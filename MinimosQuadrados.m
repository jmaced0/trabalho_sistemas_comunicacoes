MSGID = 'MATLAB:nearlySingularMatrix';
warning('off', MSGID);

load('canais_aleatorios.mat');
load('min_phase_channels.mat');

N_canal_fase = size(matrix,1);
N_canal = 1e3;
N = 1e4; %numero de amostras
M = 8; %tamanho do vetor transmitido
L = 15; %tamanho do canal
delta_min = (L-1)/2 ; %redundancia mínima
delta_max = L-1;
%delta = L-1; %redundancia maxima
faixa_SNR = 0:25;


%Transmissão

%sinal BPSK
sinal_BPSK = randi([0 1], N, 1);
sinal_BPSK(sinal_BPSK == 0) = -1;

%sinal QAM-4
sinal_QAM4 = randi([0 3], N, 1);
sinal_QAM4(sinal_QAM4 == 0) = (1+1i)*sqrt(2)/2;
sinal_QAM4(sinal_QAM4 == 1) = (1-1i)*sqrt(2)/2;
sinal_QAM4(sinal_QAM4 == 2) = (-1-1i)*sqrt(2)/2;
sinal_QAM4(sinal_QAM4 == 3) = (-1+1i)*sqrt(2)/2;

%h_comparacao = zeros(41,1);
h_comp_min = zeros(L,1);

h_comp_min(1) = 0.77+2.38i;
h_comp_min(9) = 1.58i;
h_comp_min(10) = -0.358;
h_comp_min(11) = -0.567i;
h_comp_min(14) = 0.5;
h_comp_min(15) = 0.1;

H_comp = convmtx(flip(h_comp_min),M);

%delta varia entre [7 e 14], 7 foi o caso anterior,
%então escolhemos 14 11 e 8

H_comp_delta_8 = H_comp(L-8:end,:);
H_comp_delta_8 = H_comp_delta_8*sqrt(M+8)/norm(H_comp_delta_8,'fro');

H_comp_delta_11 = H_comp(L-11:end,:);
H_comp_delta_11 = H_comp_delta_11*sqrt(M+11)/norm(H_comp_delta_11,'fro');

H_comp_delta_14 = H_comp(L-14:end,:);
H_comp_delta_14 = H_comp_delta_14*sqrt(M+14)/norm(H_comp_delta_14,'fro');

H_comp = H_comp*sqrt(22)/norm(H_comp,'fro');


%%variáveis auxiliares detecção de erros

%canal de comparacao
erros_comp_LS_delta_8_BPSK = zeros(length(faixa_SNR),1);
erros_comp_LS_delta_8_QAM4 = zeros(length(faixa_SNR),1);
erros_comp_comp_LS_delta_11_BPSK = zeros(length(faixa_SNR),1);
erros_comp_LS_delta_11_QAM4 = zeros(length(faixa_SNR),1);
erros_comp_LS_delta_14_BPSK = zeros(length(faixa_SNR),1);
erros_comp_LS_delta_14_QAM4 = zeros(length(faixa_SNR),1);


%canais aleatorios
erros_LS_delta_8_BPSK = zeros(length(faixa_SNR),1);
erros_LS_delta_8_QAM4 = zeros(length(faixa_SNR),1);
erros_LS_delta_11_BPSK = zeros(length(faixa_SNR),1);
erros_LS_delta_11_QAM4 = zeros(length(faixa_SNR),1);
erros_LS_delta_14_BPSK = zeros(length(faixa_SNR),1);
erros_LS_delta_14_QAM4 = zeros(length(faixa_SNR),1);

%canais de fase minima
erros_fase_LS_BPSK = zeros(length(faixa_SNR),1);
erros_fase_LS_QAM4 = zeros(length(faixa_SNR),1);

%sinais
u_BPSK  = reshape(sinal_BPSK, M,floor(N/M));
u_QAM4  = reshape(sinal_QAM4, M,floor(N/M));


repeticao = 5e2;
% for SNR_db = faixa_SNR
%     tic
%     SNR = db2pow(SNR_db);
%     var_ruido = 1/SNR;
%     fprintf('SNR %ddB \n',SNR_db);
% 
%     for k = 1:repeticao
%         
%         %canal especificado
%         y_BPSK = H_comp*u_BPSK + realizar_ruido(var_ruido,L+M-1,floor(N/M),1);
%         y_QAM4 = H_comp*u_QAM4 + realizar_ruido(var_ruido,L+M-1,floor(N/M),1);
%         
%         %delta 8
%         x_est_LS_QAM4 = decode(H_comp_delta_8\y_QAM4(L-8:end,:),'QAM4');
%         erros_comp_LS_delta_8_QAM4(SNR_db+1) = erros_comp_LS_delta_8_QAM4(SNR_db+1) + nnz(u_QAM4 - x_est_LS_QAM4);
% 
%         x_est_LS_BPSK = decode(H_comp_delta_8\y_BPSK(L-8:end,:),'BPSK');
%         erros_comp_LS_delta_8_BPSK(SNR_db+1) = erros_comp_LS_delta_8_BPSK(SNR_db+1) + nnz(u_BPSK - x_est_LS_BPSK);
%         
%         %delta 11
%         x_est_LS_QAM4 = decode(H_comp_delta_11\y_QAM4(L-11:end,:),'QAM4');
%         erros_comp_LS_delta_11_QAM4(SNR_db+1) = erros_comp_LS_delta_11_QAM4(SNR_db+1) + nnz(u_QAM4 - x_est_LS_QAM4);
% 
%         x_est_LS_BPSK = decode(H_comp_delta_11\y_BPSK(L-11:end,:),'BPSK');
%         erros_comp_LS_delta_11_BPSK(SNR_db+1) = erros_comp_LS_delta_11_BPSK(SNR_db+1) + nnz(u_BPSK - x_est_LS_BPSK);
%         
%         %delta 14
%         x_est_LS_QAM4 = decode(H_comp_delta_14\y_QAM4(L-14:end,:),'QAM4');
%         erros_comp_LS_delta_14_QAM4(SNR_db+1) = erros_comp_LS_delta_14_QAM4(SNR_db+1) + nnz(u_QAM4 - x_est_LS_QAM4);
% 
%         x_est_LS_BPSK = decode(H_comp_delta_14\y_BPSK(L-14:end,:),'BPSK');
%         erros_comp_LS_delta_14_BPSK(SNR_db+1) = erros_comp_LS_delta_14_BPSK(SNR_db+1) + nnz(u_BPSK - x_est_LS_BPSK);
%     end
%     toc
%      if (erros_comp_LS_delta_8_QAM4==0) %pior caso
%        break
%    end
% end
% 
% 
% erros_comp_LS_delta_8_BPSK = erros_comp_LS_delta_8_BPSK / (N*repeticao);
% erros_comp_LS_delta_8_QAM4 = erros_comp_LS_delta_8_QAM4 / (N*repeticao);
% 
% erros_comp_LS_delta_11_BPSK = erros_comp_LS_delta_11_BPSK / (N*repeticao);
% erros_comp_LS_delta_11_QAM4 = erros_comp_LS_delta_11_QAM4 / (N*repeticao);
% 
% erros_comp_LS_delta_14_BPSK = erros_comp_LS_delta_14_BPSK / (N*repeticao);
% erros_comp_LS_delta_14_QAM4 = erros_comp_LS_delta_14_QAM4 / (N*repeticao);
% 
% %plotar
% 
% %canal comparacação
% figure;
% plot(faixa_SNR,erros_comp_LS_delta_8_BPSK,'-','Color','b');
% hold on;
% plot(faixa_SNR,erros_comp_LS_delta_8_QAM4,'-','Color','r');
% hold on;
% 
% plot(faixa_SNR,erros_comp_LS_delta_11_BPSK,'--','Color','b');
% hold on;
% plot(faixa_SNR,erros_comp_LS_delta_11_QAM4,'--','Color','r');
% hold on;
% 
% plot(faixa_SNR,erros_comp_LS_delta_14_BPSK,'-.','Color','b');
% hold on;
% plot(faixa_SNR,erros_comp_LS_delta_14_QAM4,'-.','Color','r');
% hold on;
% 
% leg = legend('Ref: $\delta = 8$ \ BPSK',...
%     'Ref: $\delta = 8$ \ QAM4',...
%     'Ref: $\delta = 11$ \ BPSK',...
%     'Ref: $\delta = 11$ \ QAM4',...
%     'Ref: $\delta = 14$ \ BPSK',...
%     'Ref: $\delta = 14$ \ QAM4');
% 
% set(leg,'Interpreter', 'latex','Location','east','FontSize',14)
% 
% grid on;
% set(gca, 'YScale', 'log');
% xlabel('SNR (dB)');
% ylabel('Probabilidade de erro');







% for SNR_db = faixa_SNR
%     tic
%     SNR = db2pow(SNR_db);
%     var_ruido = 1/SNR;
%     fprintf('SNR %ddB \n',SNR_db);
%     %canais aleatorios
%     for n = 1:N_canal
%         H = squeeze(canais_aleatorios(n,:,:));
%         
%         H_delta_8 = H(L-8:end,:);
%         H_delta_8 = H_delta_8*sqrt(M+8)/norm(H_delta_8,'fro');
%         
%         H_delta_11 = H(L-11:end,:);
%         H_delta_11 = H_delta_11*sqrt(M+11)/norm(H_delta_11,'fro');
%         
%         H_delta_14 = H(L-14:end,:);
%         H_delta_14 = H_delta_14*sqrt(M+14)/norm(H_delta_14,'fro');
%                 
%         y_BPSK = H*u_BPSK + realizar_ruido(var_ruido,L+M-1,floor(N/M),1);
%         y_QAM4 = H*u_QAM4 + realizar_ruido(var_ruido,L+M-1,floor(N/M),1);
%         
%         H = H*sqrt(M+L-1)/norm(H,'fro');
%         
%         %delta 8
%         x_est_LS_QAM4 = decode(H_delta_8\y_QAM4(L-8:end,:),'QAM4');
%         erros_LS_delta_8_QAM4(SNR_db+1) = erros_LS_delta_8_QAM4(SNR_db+1) + nnz(u_QAM4 - x_est_LS_QAM4);
% 
%         x_est_LS_BPSK = decode(H_delta_8\y_BPSK(L-8:end,:),'BPSK');
%         erros_LS_delta_8_BPSK(SNR_db+1) = erros_LS_delta_8_BPSK(SNR_db+1) + nnz(u_BPSK - x_est_LS_BPSK);
%         
%         %delta 11
%         x_est_LS_QAM4 = decode(H_delta_11\y_QAM4(L-11:end,:),'QAM4');
%         erros_LS_delta_11_QAM4(SNR_db+1) = erros_LS_delta_11_QAM4(SNR_db+1) + nnz(u_QAM4 - x_est_LS_QAM4);
% 
%         x_est_LS_BPSK = decode(H_delta_11\y_BPSK(L-11:end,:),'BPSK');
%         erros_LS_delta_11_BPSK(SNR_db+1) = erros_LS_delta_11_BPSK(SNR_db+1) + nnz(u_BPSK - x_est_LS_BPSK);
%         
%         %delta 14
%         x_est_LS_QAM4 = decode(H_delta_14\y_QAM4(L-14:end,:),'QAM4');
%         erros_LS_delta_14_QAM4(SNR_db+1) = erros_LS_delta_14_QAM4(SNR_db+1) + nnz(u_QAM4 - x_est_LS_QAM4);
% 
%         x_est_LS_BPSK = decode(H_delta_14\y_BPSK(L-14:end,:),'BPSK');
%         erros_LS_delta_14_BPSK(SNR_db+1) = erros_LS_delta_14_BPSK(SNR_db+1) + nnz(u_BPSK - x_est_LS_BPSK);
%     end
%     toc
%       if (erros_LS_delta_8_QAM4==0) %pior caso
%         break
%     end
%     
% end
% erros_LS_delta_8_BPSK = erros_LS_delta_8_BPSK / (N*N_canal);
% erros_LS_delta_8_QAM4 = erros_LS_delta_8_QAM4 / (N*N_canal);
% erros_LS_delta_11_BPSK = erros_LS_delta_11_BPSK / (N*N_canal);
% erros_LS_delta_11_QAM4 = erros_LS_delta_11_QAM4 / (N*N_canal);
% erros_LS_delta_14_BPSK = erros_LS_delta_14_BPSK / (N*N_canal);
% erros_LS_delta_14_QAM4 = erros_LS_delta_14_QAM4 / (N*N_canal);
% 
% 
% %plotar
% 
% %canal comparacação
% figure;
% plot(faixa_SNR,erros_LS_delta_8_BPSK,'-','Color','b');
% hold on;
% plot(faixa_SNR,erros_LS_delta_8_QAM4,'-','Color','r');
% hold on;
% 
% plot(faixa_SNR,erros_LS_delta_11_BPSK,'--','Color','b');
% hold on;
% plot(faixa_SNR,erros_LS_delta_11_QAM4,'--','Color','r');
% hold on;
% 
% plot(faixa_SNR,erros_LS_delta_14_BPSK,'-.','Color','b');
% hold on;
% plot(faixa_SNR,erros_LS_delta_14_QAM4,'-.','Color','r');
% hold on;
% 
% leg = legend('Aleatorios: $\delta = 8$ \ BPSK',...
%     'Aleatorios: $\delta = 8$ \ QAM4',...
%     'Aleatorios: $\delta = 11$ \ BPSK',...
%     'Aleatorios: $\delta = 11$ \ QAM4',...
%     'Aleatorios: $\delta = 14$ \ BPSK',...
%     'Aleatorios: $\delta = 14$ \ QAM4');
% 
% set(leg,'Interpreter', 'latex','Location','east','FontSize',14)
% 
% grid on;
% set(gca, 'YScale', 'log');
% xlabel('SNR (dB)');
% ylabel('Probabilidade de erro');


repeticao = 20;

for SNR_db = faixa_SNR
    tic
    SNR = db2pow(SNR_db);
    var_ruido = 1/SNR;
    fprintf('SNR %ddB \n',SNR_db);
    for k=1:repeticao   
        for n = 1:N_canal_fase
            H = squeeze(matrix(n,:,:));
            
            H_delta_8 = H(L-8:end,:);
            H_delta_8 = H_delta_8*sqrt(M+8)/norm(H_delta_8,'fro');
            
            H_delta_11 = H(L-11:end,:);
            H_delta_11 = H_delta_11*sqrt(M+11)/norm(H_delta_11,'fro');
            
            H_delta_14 = H(L-14:end,:);
            H_delta_14 = H_delta_14*sqrt(M+14)/norm(H_delta_14,'fro');
            
            y_BPSK = H*u_BPSK + realizar_ruido(var_ruido,L+M-1,floor(N/M),1);
            y_QAM4 = H*u_QAM4 + realizar_ruido(var_ruido,L+M-1,floor(N/M),1);
            
            H = H*sqrt(M+L-1)/norm(H,'fro');
            
            %delta 8
            x_est_LS_QAM4 = decode(H_delta_8\y_QAM4(L-8:end,:),'QAM4');
            erros_LS_delta_8_QAM4(SNR_db+1) = erros_LS_delta_8_QAM4(SNR_db+1) + nnz(u_QAM4 - x_est_LS_QAM4);
            
            x_est_LS_BPSK = decode(H_delta_8\y_BPSK(L-8:end,:),'BPSK');
            erros_LS_delta_8_BPSK(SNR_db+1) = erros_LS_delta_8_BPSK(SNR_db+1) + nnz(u_BPSK - x_est_LS_BPSK);
            
            %delta 11
            x_est_LS_QAM4 = decode(H_delta_11\y_QAM4(L-11:end,:),'QAM4');
            erros_LS_delta_11_QAM4(SNR_db+1) = erros_LS_delta_11_QAM4(SNR_db+1) + nnz(u_QAM4 - x_est_LS_QAM4);
            
            x_est_LS_BPSK = decode(H_delta_11\y_BPSK(L-11:end,:),'BPSK');
            erros_LS_delta_11_BPSK(SNR_db+1) = erros_LS_delta_11_BPSK(SNR_db+1) + nnz(u_BPSK - x_est_LS_BPSK);
            
            %delta 14
            x_est_LS_QAM4 = decode(H_delta_14\y_QAM4(L-14:end,:),'QAM4');
            erros_LS_delta_14_QAM4(SNR_db+1) = erros_LS_delta_14_QAM4(SNR_db+1) + nnz(u_QAM4 - x_est_LS_QAM4);
            
            x_est_LS_BPSK = decode(H_delta_14\y_BPSK(L-14:end,:),'BPSK');
            erros_LS_delta_14_BPSK(SNR_db+1) = erros_LS_delta_14_BPSK(SNR_db+1) + nnz(u_BPSK - x_est_LS_BPSK);
        end
    end
    toc
    
    if (erros_LS_delta_8_QAM4(SNR_db+1)==0) %pior caso
        break
    end
    
end
erros_LS_delta_8_BPSK = erros_LS_delta_8_BPSK / (N*N_canal_fase*repeticao);
erros_LS_delta_8_QAM4 = erros_LS_delta_8_QAM4 / (N*N_canal_fase*repeticao);
erros_LS_delta_11_BPSK = erros_LS_delta_11_BPSK / (N*N_canal_fase*repeticao);
erros_LS_delta_11_QAM4 = erros_LS_delta_11_QAM4 / (N*N_canal_fase*repeticao);
erros_LS_delta_14_BPSK = erros_LS_delta_14_BPSK / (N*N_canal_fase*repeticao);
erros_LS_delta_14_QAM4 = erros_LS_delta_14_QAM4 / (N*N_canal_fase*repeticao);


%plotar

%canal comparacação
figure;
plot(faixa_SNR,erros_LS_delta_8_BPSK,'-','Color','b');
hold on;
plot(faixa_SNR,erros_LS_delta_8_QAM4,'-','Color','r');
hold on;

plot(faixa_SNR,erros_LS_delta_11_BPSK,'--','Color','b');
hold on;
plot(faixa_SNR,erros_LS_delta_11_QAM4,'--','Color','r');
hold on;

plot(faixa_SNR,erros_LS_delta_14_BPSK,'-.','Color','b');
hold on;
plot(faixa_SNR,erros_LS_delta_14_QAM4,'-.','Color','r');
hold on;

leg = legend('fase minima: $\delta = 8$ \ BPSK',...
    'fase minima: $\delta = 8$ \ QAM4',...
    'fase minima: $\delta = 11$ \ BPSK',...
    'fase minima: $\delta = 11$ \ QAM4',...
    'fase minima: $\delta = 14$ \ BPSK',...
    'fase minima: $\delta = 14$ \ QAM4');

set(leg,'Interpreter', 'latex','Location','east','FontSize',14)

grid on;
set(gca, 'YScale', 'log');
xlabel('SNR (dB)');
ylabel('Probabilidade de erro');
