clear all;
close all;
N = 256;
x = randn( 1, N ); %wol� mie� w wierszu jak na wyk�adzie
X1 = fft(x); % oryginalne DFT
X2 = dit(x); % DFT ,,sklejane'' z dw�ch po��wek
X3 = recDit(x); % rekursywne FFT


mean( abs(X1-X2) ) % b��d
mean( abs(X1-X3) ) % b��d

