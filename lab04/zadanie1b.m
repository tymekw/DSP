clear all;
close all;
N = 256;
x = randn( 1, N ); %wolê mieæ w wierszu jak na wyk³adzie
X1 = fft(x); % oryginalne DFT
X2 = dit(x); % DFT ,,sklejane'' z dwóch po³ówek
X3 = recDit(x); % rekursywne FFT


mean( abs(X1-X2) ) % b³¹d
mean( abs(X1-X3) ) % b³¹d

