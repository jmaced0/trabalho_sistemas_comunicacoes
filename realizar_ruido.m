function Z = realizar_ruido(sigma,i, j, k)
%sigma = desvio padrao do ruido
%Z = vetor contendo (N,M) realizacoes do ruido
%N = quantidade de realizacoes
%M = tamanho do pulso

%definicao do ruido branco
Z = normrnd(0,sqrt(sigma/2),[i,j,k]) + 1i*normrnd(0,sqrt(sigma/2),[i,j,k]);


