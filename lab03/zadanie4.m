clc;
clear all;

load("lab_03.mat");
pytajniki = mod(305525, 16)+1;
%x_6
K = 8;
N=512;
M=32;

%na ramki
for m = 0:K-1
    ramki(m+1,:)=x_6(m*(N+M)+1+M:m*(N+M)+M+N);
end

%N punktowe DFT (FFT) dla kazdej ramki
for m = 0:K-1
   dft(m+1,:) = fft(ramki(m+1,:))./N;
end

f=1:512;
%rysunki
for m = 0:K-1
   figure;
   plot(f,20*log10(abs(dft(m+1,:))));
end

%fs = ?
%f0 = 1/T = 1/(N/fs)
%f k-tej probki => k*f0


