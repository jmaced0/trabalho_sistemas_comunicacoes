function out = decode(data, codification)
    out = zeros(size(data));

    if (strcmp (codification, 'BPSK'))
        
        out(real(data) >= 0) = 1;
        out(real(data) < 0) = -1;
        
    elseif (strcmp (codification, 'QAM4'))
        
        out((real(data) >= 0) & (imag(data) >= 0)) = (1+1i)*sqrt(2)/2;
        out((real(data) >= 0) & (imag(data) < 0 )) = (1-1i)*sqrt(2)/2;
        out((real(data) < 0)  & (imag(data) < 0 )) = (-1-1i)*sqrt(2)/2;
        out((real(data) < 0)  & (imag(data) >= 0)) = (-1+1i)*sqrt(2)/2;
        
    else
        disp('codification is not available')
    end
end


