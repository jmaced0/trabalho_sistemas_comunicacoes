function [ canal ] = gerar_canal(L)
    canal = (1/sqrt(2))*(randn(L,1) +1i*randn(L,1));
end