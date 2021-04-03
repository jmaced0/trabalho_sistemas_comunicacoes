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
delta = L-2;
P = M+2*delta-L+1; %largura da matriz H
%delta = L-1; %redundancia maxima
faixa_SNR = 0:3;


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
H_comp_min = H_comp;
H_comp = H_comp*sqrt(L+M-1)/norm(H_comp,'fro');

%%variáveis auxiliares detecção de erros

%canal de comparacao
erros_comp_min_BPSK = zeros(length(faixa_SNR),1);
erros_comp_min_QAM4 = zeros(length(faixa_SNR),1);
%canais aleatorios
erros_min_BPSK = zeros(length(faixa_SNR),1);
erros_min_QAM4 = zeros(length(faixa_SNR),1);
%canais de fase minima
erros_fase_min_BPSK = zeros(length(faixa_SNR),1);
erros_fase_min_QAM4 = zeros(length(faixa_SNR),1);

%sinais
u_BPSK  = reshape(sinal_BPSK, M,floor(N/M));
u_QAM4  = reshape(sinal_QAM4, M,floor(N/M));




% erros_comp = zeros(15,1);
% erros_random = zeros(15,1);
%
%
% for i=1:L
%     H_comp_min_loop = H_comp(i:i+M-1,:);
%     H_comp_min_loop = H_comp_min_loop*sqrt(M)/norm(H_comp_min_loop,'fro');
%
%     y_QAM4 = H_comp*u_QAM4+ realizar_ruido(var_ruido,L+M-1,floor(N/M),1);
%     x_est_QAM4_min = decode(H_comp_min_loop\y_QAM4(end-7:end,:),'QAM4');
%     erros_comp(i) =  erros_comp(i)  + nnz(u_QAM4 - x_est_QAM4_min);
%
%     for k = 1:N_canal
%         H = squeeze(canais_aleatorios(k,:,:));
%         H_min = H(i:i+M-1,:);
%         H = H*sqrt(M+L-1)/norm(H,'fro');
%         H_min = H_min*sqrt(M)/norm(H_min,'fro');
%
%         y_QAM4 = H*u_QAM4+ realizar_ruido(var_ruido,L+M-1,floor(N/M),1);
%         x_est_QAM4_min = decode(H_min\y_QAM4(end-7:end,:),'QAM4');
%         erros_random(i) = erros_random(i) + nnz(u_QAM4 - x_est_QAM4_min);
%     end
%     erros_random(i) = erros_random(i) / N_canal;
% end
%
% [ max_comp, max_index] = min(erros_comp);
% best_slice_comp = max_index;
%
% [max_random , max_index] = min(erros_random);
% best_slice_random = max_index;

%apenas para não precisar rodar o codigo acima novamente

H_comp_min = H_comp_min(L:L+M-1,:);
H_comp_min = H_comp_min*sqrt(M)/norm(H_comp_min,'fro');

% repeticao = 5e2;
% 
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
%         
%         
%         %inversao de H0
%         x_est_QAM4_min = decode(H_comp_min\y_QAM4(end-7:end,:),'QAM4');
%         erros_comp_min_QAM4(SNR_db+1) =  erros_comp_min_QAM4(SNR_db+1) + nnz(u_QAM4 - x_est_QAM4_min);
%         
%         x_est_BPSK_min = decode(H_comp_min\y_BPSK(end-7:end,:),'BPSK');
%         erros_comp_min_BPSK(SNR_db+1) =  erros_comp_min_BPSK(SNR_db+1) + nnz(u_BPSK - x_est_BPSK_min);
%         
%         
%         
%     end
%     toc
% end


% erros_comp_min_BPSK = erros_comp_min_BPSK / (N*repeticao);
% erros_comp_min_QAM4 = erros_comp_min_QAM4 / (N*repeticao);
% 
% for SNR_db = faixa_SNR
%     tic
%     SNR = db2pow(SNR_db);
%     var_ruido = 1/SNR;
%     fprintf('SNR %ddB \n',SNR_db);
%     %canais aleatorios
%     for n = 1:N_canal
%         H = squeeze(canais_aleatorios(n,:,:));
%         H_min = H(L:L+M-1,:);
%         H = H*sqrt(M+L-1)/norm(H,'fro');
%         H_min = H_min*sqrt(M)/norm(H_min,'fro');
%         
%         
%         y_BPSK = H*u_BPSK + realizar_ruido(var_ruido,L+M-1,floor(N/M),1);
%         y_QAM4 = H*u_QAM4 + realizar_ruido(var_ruido,L+M-1,floor(N/M),1);
%         
%         
%         
%         %inversao de H0
%         x_est_QAM4_min = decode(H_min\y_QAM4(end-7:end,:),'QAM4');
%         erros_min_QAM4(SNR_db+1) = erros_min_QAM4(SNR_db+1) + nnz(u_QAM4 - x_est_QAM4_min);
%         
%         x_est_BPSK_min = decode(H_min\y_BPSK(end-7:end,:),'BPSK');
%         erros_min_BPSK(SNR_db+1) = erros_min_BPSK(SNR_db+1) + nnz(u_BPSK - x_est_BPSK_min);
%         
%         
%     end
%     toc
%     
% end
% erros_min_BPSK = erros_min_BPSK / (N*N_canal);
% erros_min_QAM4 = erros_min_QAM4 / (N*N_canal);


repeticao = 1e3;%1e2;
for SNR_db = faixa_SNR
    tic
    SNR = db2pow(SNR_db);
    var_ruido = 1/SNR;
    
    for k = 1:repeticao
        %canais fase minima
        for n = 1:1%N_canal_fase
            H = squeeze(matrix(n,:,:));
            H_min = H(L:L+M-1,:);
            H = H*sqrt(M+L-1)/norm(H,'fro');
            H_min = H_min*sqrt(M)/norm(H_min,'fro');
            
            
            y_BPSK = H*u_BPSK + realizar_ruido(var_ruido,L+M-1,floor(N/M),1);
            y_QAM4 = H*u_QAM4 + realizar_ruido(var_ruido,L+M-1,floor(N/M),1);
           
            
            %inversao de H0
            x_est_QAM4_min = decode(H_min\y_QAM4(end-7:end,:),'QAM4');
            erros_fase_min_QAM4(SNR_db+1) = erros_fase_min_QAM4(SNR_db+1) + nnz(u_QAM4 - x_est_QAM4_min);
            
            x_est_BPSK_min = decode(H_min\y_BPSK(end-7:end,:),'BPSK');
            erros_fase_min_BPSK(SNR_db+1) = erros_fase_min_BPSK(SNR_db+1) + nnz(u_BPSK - x_est_BPSK_min);
            
            
        end
        
    end
    fprintf('SNR %ddB \n',SNR_db);
    toc
end
%calcular média dos erros


erros_fase_min_BPSK = erros_min_BPSK / (N*repeticao);%(N*N_canal_fase*repeticao);
erros_fase_min_QAM4 = erros_min_QAM4 / (N*repeticao);%(N*N_canal_fase*repeticao);

%plotar

%canal comparacação
plot(faixa_SNR,erros_comp_min_BPSK,'-','Color','b');
hold on;
plot(faixa_SNR,erros_comp_min_QAM4,'-','Color','r');
hold on;

%canais aleatórios
plot(faixa_SNR,erros_min_BPSK,'--','Color','b');
hold on;
plot(faixa_SNR,erros_min_QAM4,'--','Color','r');
hold on;

%canais fase minima
plot(faixa_SNR,erros_fase_min_BPSK,'-.','Color','b');
hold on;
plot(faixa_SNR,erros_fase_min_QAM4,'-.','Color','r');
hold on;

leg = legend('Canal especificado: BPSK',...
    'Canal especificado: QAM4',...
    'Canais aleatórios: BPSK',...
    'Canais aleatórios: QAM4',...
    'Canais fase mínima: BPSK',...
    'Canais fase mínima: QAM4');

set(leg,'FontSize',14,'Location','east' )

grid on;
set(gca, 'YScale', 'log');
xlabel('SNR (dB)');
ylabel('Probabilidade de erro');




