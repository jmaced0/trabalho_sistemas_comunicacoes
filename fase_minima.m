MSGID = 'MATLAB:nearlySingularMatrix';
warning('off', MSGID);

load('canais_aleatorios.mat');

N_canal = 1e4;
N = 1e3; %numero de amostras
M = 8; %tamanho do vetor transmitido
L = 15; %tamanho do canal
delta_min = (L-1)/2 ; %redundancia mínima
delta = L-2;
P = M+2*delta-L+1; %largura da matriz H
%delta = L-1; %redundancia maxima
faixa_SNR = 0:25;


                        %Transmissão

%sinal BPSK
sinal_BPSK = randi([0 1], N, 1);
sinal_BPSK(sinal_BPSK == 0) = -1;



h=zeros(L,1);
h(1) = 0.77+2.38i;
h(9) = 1.58i;
h(10) = -0.358;
h(11) = -0.567i;
h(14) = 0.5;
h(15) = 0.1;

ts = 1e-4;

min_phase_channels = [];

for n = 1:N_canal
    aux = flip(canais_aleatorios(n,1:L,1));
    filter_order = find(aux,1,'last') ;
    sys = tf(aux,1, ts, 'Variable','z^-1');
    [p, z] = pzmap(sys);
    if sum(abs(z) < 1) == filter_order-1
        min_phase_channels  = [min_phase_channels; aux];
    end
end

[m,n] = size(min_phase_channels);
matrix = zeros(m,L+M-1,M);
for k = 1:m
    matrix(k,:,:) = convmtx(flip(min_phase_channels(k,:)'),M);
end

save('min_phase_channels.mat','matrix');

ts = 1e-5;
sys_comp_min = tf(h',1, ts, 'Variable','z^-1');
[p_comp, z_comp] = pzmap(sys_comp_min);

aux = min_phase_channels(2,:);
sys = tf(aux,1, ts, 'Variable','z^-1');
[p_min_phase, z_min_phase] = pzmap(sys);

aux = canais_aleatorios(1,1:L,1);
sys = tf(aux,1, ts, 'Variable','z^-1');
[p, z] = pzmap(sys);

figure;

subplot(1,3,1)
zplane(z_comp, p_comp);
leg = legend('canal especificado');
set(leg,'FontSize',14)

grid
ylim([-1.05, 1.25]);
xlim([-1.05, 1.05]);

subplot(1,3,2)
zplane(z, p);
leg = legend('canal aleatório');
set(leg,'FontSize',14)
grid
ylim([-1.05, 1.25]);
xlim([-1.05, 1.05]);

subplot(1,3,3)
zplane( z_min_phase,p_min_phase);
leg = legend('canal de fase mínima');
set(leg,'FontSize',14)
grid
ylim([-1.05, 1.25]);
xlim([-1.05, 1.05]);

save('min_phase_channels.mat','matrix');
