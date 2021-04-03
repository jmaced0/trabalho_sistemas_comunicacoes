MSGID = 'MATLAB:nearlySingularMatrix';
warning('off', MSGID);
MSGID = 'MATLAB:singularMatrix';
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
h_ref_min = zeros(L,1);

h_ref_min(1) = 0.77+2.38i;
h_ref_min(9) = 1.58i;
h_ref_min(10) = -0.358;
h_ref_min(11) = -0.567i;
h_ref_min(14) = 0.5;
h_ref_min(15) = 0.1;
h_ref_min = h_ref_min/norm(h_ref_min);

H_ref = convmtx(h_ref_min,L);
H_ref = H_ref*sqrt(M+L-1)/norm(H_ref,'fro');

%%variáveis auxiliares detecção de erros

%canal de comparacao
erros_ref_BPSK = zeros(length(faixa_SNR),1);
erros_ref_QAM4 = zeros(length(faixa_SNR),1);
%canais aleatorios
erros_rnd_BPSK = zeros(length(faixa_SNR),1);
erros_rnd_QAM4 = zeros(length(faixa_SNR),1);
%canais de fase minima
erros_fase_BPSK = zeros(length(faixa_SNR),1);
erros_fase_QAM4 = zeros(length(faixa_SNR),1);

%sinais
s_BPSK  = reshape(sinal_BPSK, M,floor(N/M));
s_QAM4  = reshape(sinal_QAM4, M,floor(N/M));


F = dftmtx(M+L-1);
lambda = sqrt(M) *fft(h_ref_min(1:M));
Lambda = diag(lambda);
% C = F'*Lambda*F;
% C(abs(C) < 1e-7)=0;
% C = C*sqrt(22)/norm(C,'fro');

%C2 = bsxfun(@circshift, h_ref_min(1:M),0:M-1 );


%calcular IDFT do sinal
s_barra_BPSK = ifft(s_BPSK,[],1);%ifft divide a potencia do sinal pelo tamnho do vetor
s_barra_QAM4 = ifft(s_QAM4,[],1);%dividir por 8 para ajustar a SNR

%adicionar prefico ciclico
u_BPSK = [s_barra_BPSK; s_barra_BPSK(1:M-1,:) ];
u_QAM4 = [s_barra_QAM4; s_barra_QAM4(1:M-1,:) ];

Psinal = mean(var(u_QAM4)); 

% repeticao = 1e3;
% 
% for SNR_db = faixa_SNR
%     tic
%     SNR = db2pow(SNR_db);
%     var_ruido = Psinal/SNR;
%     fprintf('SNR %ddB \n',SNR_db);
%     
%     for k = 1:repeticao
%     
%     %passar pelo canal
%     y_BPSK = H_ref*u_BPSK + realizar_ruido(var_ruido,size(H_ref,1),floor(N/M),1);
%     y_QAM4 = H_ref*u_QAM4 + realizar_ruido(var_ruido,size(H_ref,1),floor(N/M),1);
%     
%     %remover prefixo ciclico
%     y_barra_BPSK = y_BPSK(1:M,:);
%     y_barra_QAM4 = y_QAM4(1:M,:);
%     
%     %passar pela dft
%     y_barra_BPSK = fft(y_barra_BPSK,[],1);
%     y_barra_QAM4 = fft(y_barra_QAM4,[],1);
%     
%     %detecção
%     s_estimado_BPSK = Lambda\y_barra_BPSK;
%     s_estimado_QAM4 = Lambda\y_barra_QAM4;    
%     
%     erros_ref_BPSK(SNR_db+1) =  erros_ref_BPSK(SNR_db+1) + nnz(decode(s_estimado_BPSK,'BPSK') - s_BPSK);
%     erros_ref_QAM4(SNR_db+1) =  erros_ref_QAM4(SNR_db+1) + nnz(decode(s_estimado_QAM4,'QAM4') - s_QAM4);
%      
%     end
%     erros_ref_BPSK(SNR_db+1) = erros_ref_BPSK(SNR_db+1)/(N*repeticao);
%     erros_ref_QAM4(SNR_db+1) = erros_ref_QAM4(SNR_db+1)/(N*repeticao);
%     
%     if erros_ref_QAM4(SNR_db+1) == 0
%         break
%     end
%     
% %     erros_ref_BPSK(SNR_db+1) =  nnz(decode(s_estimado_BPSK,'BPSK') - s_BPSK)/N;
% %     erros_ref_QAM4(SNR_db+1) =  nnz(decode(s_estimado_QAM4,'QAM4') - s_QAM4)/N;
%     
%     toc
% end
% 
% figure;
% plot(faixa_SNR,erros_ref_BPSK,'-','Color','b');
% hold on;
% plot(faixa_SNR,erros_ref_QAM4,'-','Color','r');
% hold on;
% 
% leg = legend('Ref: OFDM BPSK',...
%             'Ref: OFDM QAM4');
% 
% set(leg,'Interpreter', 'latex','Location','east','FontSize',14)
% grid on;
% set(gca, 'YScale', 'log');
% xlabel('SNR (dB)');
% ylabel('Probabilidade de erro');


repeticao = 10;

for SNR_db = faixa_SNR
    tic
    SNR = db2pow(SNR_db);
    var_ruido = Psinal/SNR;
    fprintf('SNR %ddB \n',SNR_db);
    for n = 1:N_canal_fase
        h = flip(squeeze(matrix(n,1:L,1)))';
        h = h/norm(h);
        
        
        H = convmtx(h,L);
        H = H*sqrt(M+L-1)/norm(H,'fro');
        
        lambda = sqrt(M) *fft(h(1:M));
        Lambda = diag(lambda);
        
        for k = 1:repeticao
            
            %passar pelo canal
            y_BPSK = H*u_BPSK + realizar_ruido(var_ruido,size(H,1),floor(N/M),1);
            y_QAM4 = H*u_QAM4 + realizar_ruido(var_ruido,size(H,1),floor(N/M),1);
            
            %remover prefixo ciclico
            y_barra_BPSK = y_BPSK(1:M,:);
            y_barra_QAM4 = y_QAM4(1:M,:);
            
            %passar pela dft
            y_barra_BPSK = fft(y_barra_BPSK,[],1);
            y_barra_QAM4 = fft(y_barra_QAM4,[],1);
            
            %detecção
            s_estimado_BPSK = Lambda\y_barra_BPSK;
            s_estimado_QAM4 = Lambda\y_barra_QAM4;
            
            erros_fase_BPSK(SNR_db+1) =  erros_fase_BPSK(SNR_db+1) + nnz(decode(s_estimado_BPSK,'BPSK') - s_BPSK);
            erros_fase_QAM4(SNR_db+1) =  erros_fase_QAM4(SNR_db+1) + nnz(decode(s_estimado_QAM4,'QAM4') - s_QAM4);
            
        end
        erros_fase_BPSK(SNR_db+1) = erros_fase_BPSK(SNR_db+1)/(N*repeticao*N_canal_fase);
        erros_fase_QAM4(SNR_db+1) = erros_fase_QAM4(SNR_db+1)/(N*repeticao*N_canal_fase);
        
        
        
    end
    
    if (erros_fase_QAM4(SNR_db+1) < 1e-7)
        break
    end
    toc
end

figure;
plot(faixa_SNR,erros_fase_BPSK,'-','Color','b');
hold on;
plot(faixa_SNR,erros_fase_QAM4,'-','Color','r');
hold on;

leg = legend('fase minima: OFDM BPSK',...
            'fase minima: OFDM QAM4');

set(leg,'Interpreter', 'latex','Location','northeast','FontSize',14)
grid on;
set(gca, 'YScale', 'log');
xlabel('SNR (dB)');
ylabel('Probabilidade de erro');




% repeticao = 2;
% 
% for SNR_db = faixa_SNR
%     tic
%     SNR = db2pow(SNR_db);
%     var_ruido = Psinal/SNR;
%     fprintf('SNR %ddB \n',SNR_db);
%     for n = 1:N_canal
%         h = flip(squeeze(canais_aleatorios(n,1:L,1)))';
%         h = h/norm(h);
%         
%         
%         H = convmtx(h,L);
%         H = H*sqrt(M+L-1)/norm(H,'fro');
%         
%         lambda = sqrt(M) *fft(h(1:M));
%         Lambda = diag(lambda);
%         
%         for k = 1:repeticao
%             
%             %passar pelo canal
%             y_BPSK = H*u_BPSK + realizar_ruido(var_ruido,size(H,1),floor(N/M),1);
%             y_QAM4 = H*u_QAM4 + realizar_ruido(var_ruido,size(H,1),floor(N/M),1);
%             
%             %remover prefixo ciclico
%             y_barra_BPSK = y_BPSK(1:M,:);
%             y_barra_QAM4 = y_QAM4(1:M,:);
%             
%             %passar pela dft
%             y_barra_BPSK = fft(y_barra_BPSK,[],1);
%             y_barra_QAM4 = fft(y_barra_QAM4,[],1);
%             
%             %detecção
%             s_estimado_BPSK = Lambda\y_barra_BPSK;
%             s_estimado_QAM4 = Lambda\y_barra_QAM4;
%             
%             erros_rnd_BPSK(SNR_db+1) =  erros_rnd_BPSK(SNR_db+1) + nnz(decode(s_estimado_BPSK,'BPSK') - s_BPSK);
%             erros_rnd_QAM4(SNR_db+1) =  erros_rnd_QAM4(SNR_db+1) + nnz(decode(s_estimado_QAM4,'QAM4') - s_QAM4);
%             
%         end
%         erros_rnd_BPSK(SNR_db+1) = erros_rnd_BPSK(SNR_db+1)/(N*repeticao*N_canal);
%         erros_rnd_QAM4(SNR_db+1) = erros_rnd_QAM4(SNR_db+1)/(N*repeticao*N_canal);
%         
%         
%         
%     end
%     if (erros_rnd_QAM4(SNR_db+1) < 1e-7)
%              break
%     end
%     toc
% end
% 
% figure;
% plot(faixa_SNR,erros_rnd_BPSK,'-','Color','b');
% hold on;
% plot(faixa_SNR,erros_rnd_QAM4,'-','Color','r');
% hold on;
% 
% leg = legend('Aleatorios: OFDM BPSK',...
%             'Aleatorios: OFDM QAM4');
% 
% set(leg,'Interpreter', 'latex','Location','northeast','FontSize',14)
% grid on;
% set(gca, 'YScale', 'log');
% xlabel('SNR (dB)');
% ylabel('Probabilidade de erro');
