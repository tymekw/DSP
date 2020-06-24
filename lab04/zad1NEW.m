clear all;
close all;
N = 256;
x = randn( N, 1 ); %w kolumnie
X1 = fft(x); % oryginalne DFT
X2 = dit1(x); % DFT ,,sklejane'' z dw�ch po��wek
X3 = recDit1(x); % rekursywne FFT
X3=X3(:);


mean( abs(X1-X2) ) % b��d
mean( abs(X1-X3) ) % b��d


fs = 256;
f = fs*(0:N-1)/N;
stem(f,abs(X3));
