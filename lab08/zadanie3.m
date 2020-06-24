clc;
clear all; close all;
T=1;
fs1 = 8000;
Nx1 = T*fs1;
t1 = (0:1:Nx1-1)*1/fs1;
f1 = 1001.2;
x1 = sin(2*pi*t1*f1);

fs2 = 32000;
Nx2 = T*fs2;
t2 = (0:1:Nx2-1)*1/fs2;
f2 = 303.1;
x2 = sin(2*pi*t2*f2);

fs3 = 48000;
Nx3 = T*fs3;
t3 = (0:1:Nx3-1)*1/fs3;
f3 = 2110.4;
x3 = sin(2*pi*t3*f3);

figure;
plot(t1,x1);
hold on;
plot(t2,x2);
plot(t3,x3);
xlim([0,5]*10^(-4));
legend("1001.2","303.1","2110.4")


dt4 = 1/fs3;
tw = (0:dt4:1-dt4);
x4w = sin(2*pi*f1*tw)+sin(2*pi*f2*tw)+sin(2*pi*f3*tw);

%upsampling 8->48; 6razy
K = length(x3)/length(x1);
xz1 = zeros(1,K*Nx1);
xz1(1:K:end) = x1;
h = K*fir1(101,1/K);

%{
fs=48000;
f=0:1:fs;
figure;
H = freqz(h,1,f,fs);
plot(f,20*log10(abs(H)));
%}



x1i = filter(h,1,xz1);

x1im = resample(x1,Nx3,Nx1);
x1ii = interp(x1,K);

x1u = x1i;  %wybranie upsampling

%porownanie mojej wersji z matlabem
figure;
plot(x1i);
hold on;
plot(x1im);
plot(x1ii);
legend("rêcznie","matlab resample","matlab interp");
title("porównanie funkcji up-sample");
xlim([0,200]);



%upsampling 32->48 ale K=1.5, slabo
%upsampling 32->96 K=3; ->
%downsampling 96 -> 48; L=2;
K = 3;
xz2 = zeros(1,K*Nx2);
xz2(1:K:end) = x2;
h = K*fir1(101,1/K); %LP filter 1/K = fs2/fs3 przepuszcza do 1/K
%x2i = filter(h,1,xz2);

x2i = conv(xz2,h);
x2i = x2i(101+1:end);
x2u = x2i;


%down
L = 2;
g = fir1(101,1/L-0.1*(1-L)); %LP filter zostawia 1/L-t¹ czest sygnalu nie przepuszca od 1/L
xd = conv(x2u,g);
x2d = xd(101+1:end);
x2ud = x2d(1:L:end);

x4 = x1u+x2ud+x3;

figure;
plot(x4w);
hold on;
plot(x4);
legend("wygenerowany","moj");

%soundsc(x4,fs3); %stuków klików nie ma
%pause;

%miksowanie rzeczywistych
clear all;
[x1,fs1] = audioread("x1.wav");
[x2,fs2] = audioread("x2.wav");
fpr = 48000;
K1 = fpr/fs1; %s³abo bo nieca³kowite
K2 = fpr/fs2;

%upsampling K2 signal x2.wav
Nx2 = length(x2);
xz2 = zeros(1,K2*Nx2);
xz2(1:K2:end) = x2;
h = K2*fir1(101,1/K2);

x2i = filter(h,1,xz2);
x2im = resample(x2,fpr,fs2);
x2ii = interp(x2,K2);
x2u = x2i;  %wybranie upsampling



%upsampling 3 times K1 then downsampling it 2 times
x1 = x1(:,2);
K = 3;
Nx1 = length(x1);
xz1 = zeros(1,K*Nx1);
xz1(1:K:end) = x1;
h = K*fir1(101,1/K);
x1i = filter(h,1,xz1);
x1u = x1i;

%down
L = 2;
g = fir1(101,1/L-0.1*(1-L));
xd = conv(x1u,g);
x1d = xd(101+1:end);
x1ud = x1d(1:L:end);

x2uz = zeros(1,length(x1ud));
x2uz(1:length(x2u)) = x2u;

%xsum = x1ud(1:length(x2u))+x2u;
xsum = x1ud + x2uz;

%soundsc(xsum,48000); % ..be happy





















