function xa=myHilbert(x)
    N = length(x); % input signal length
    X = fft(x); % signal FFT, then its modification:
    X(1)=0; %X(N/2+1)=0; % # 0 for 0 Hz and fs/2
    X(2:N/2) = -j*X(2:N/2); % # (-j) for positive frequencies
    X(N/2+2:N) = j*X(N/2+2:N); % # (+j) for negative frequencies
    xi = ifft(X); % inverse FFT of the modified spectrum
    xa = x+j*xi;
end