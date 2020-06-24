close all;
clear all;

[s2, fs] = audioread('s2.wav');
figure(1);
spectrogram(s2, 4096, 4096-512, [0:5:2000], fs); % 0,6,5,1,2

load('butter.mat');
f = 0:1:fs/2-1;

[zz,pp,kk] = bilinear(z,p,k,fs);
%bC = kk* poly(zz);
%aC = poly(pp);
[bC,aC] = zp2tf(zz,pp,kk);
HC = freqz(bC,aC,f,fs);
%plot(f,20*log10(abs(HC)),'r');
s2afterfilter = filter(bC,aC,s2);
figure(2)
spectrogram(s2afterfilter, 4096, 4096-512, [0:5:2000], fs);

N = length(s2);
T=6;
fpr = N/T;
dt = 1/fpr;
t = dt*(0:N-1);
%t = (1/N):(6/N):6;
figure(3)
plot(t,s2,'x');
hold on;
plot(t, s2afterfilter, 'o');
