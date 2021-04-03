L = 15;
M = 8;
N = 1e4;

canais_aleatorios = zeros(N,L+M-1,M);

for n = 1:N
    aux = gerar_canal(L);
    h(1) = aux(1); % primeiro tap
    h(2:L) = randi([0 1], L-1,1).* aux(2:L); %gerar matriz esparsa
    canais_aleatorios(n,:,:)  = convmtx(flip(h),M);  
end

save('canais_aleatorios.mat','canais_aleatorios');
