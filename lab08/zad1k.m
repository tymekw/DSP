clear all; close all; clc;
load('lab08_am.mat');
x = s2;
    
T = 1;
fs = 1000;
dt = 1/fs;
t = (0:fs*T -1)*dt;
A1=0.36;
A2=0.1;
A3=0.26;
%PAMI?TA? O JEDYNCE xd
x2 = 1+ A1*cos(2*pi*8*t)+A2*cos(2*pi*30*t)+A3*cos(2*pi*60*t);

xa=myHilbert(x);
AM = abs(xa);
err = abs(x2-AM);



N = length(AM);
amFourier = fft(AM)./length(AM);
amFourier(1) = 0; %amFourier(N/2+1) = 0;
f = (0:length(amFourier)-1)*fs/length(amFourier);

figure;
plot(f,abs(amFourier))
title('AM fourier transform')

figure;
plot(t,AM); title('AM(t)'); grid;

figure;
plot(t,x2);
hold on;
plot(t,err);
legend("reconstructed","error")

