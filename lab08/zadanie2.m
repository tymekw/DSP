clc;
clear all; close all;

%filtr hilberta
  M=2500; n =-M:M; MM = 2*M+1;
  hH = (1-cos(pi*n))./(pi*n); hH(M+1)=0;
  h=hH;
%okno
  %w=blackman(2*M+1)';
  w=kaiser(2*M+1,10)';
  h=h.*w;
  
%{
fs=400e3;
f=-fs/2:1:fs/2;
H = freqz(h,1,f,fs);
figure;
hold on;
plot(f,20*log10(abs(H)));
%}
%xr = x(M+1:Nx-M);
%xi = y(2*M+1:Nx);

fs = 400e3;
fc1 = 100e3;
[x1,fsx] = audioread("mowa8000.wav");
x1=x1';
fc2 = 110e3;
dA=0.25;

s=1;
for n=length(x1):-1:1
 x2(s) = x1(n);
 s=s+1;
end
clear s;


%nadpróbkowanie x1 i x2 fs/fsx razy
K = fs/fsx;

%ksi¹¿ka
if(1)
  x2i = interp(x2,K);
  x1i = interp(x1,K);
end

if(0)
%na piechotê ale cos nie dziala niestety 
  xz = zeros(1,K*length(x1));
  xz(1:K:end) = x1;
  h = K*fir1(101,1/K);
  x1 = filter(h,1,xz);
  x1i = x1;
  
  xz = zeros(1,K*length(x2));
  xz(1:K:end) = x2;
  x2 = filter(h,1,xz);
  x2i = x2;
end


%filtracja
Nx = length(x1i);
xh1 = conv(x1i,h);
%xr1 = x1(M+1:Nx-M);
xi1 = xh1(2*M+1:Nx);
%xa1 = xr1 + j*xi1;
%xh1 = abs(xa1);

xh2 = conv(x2i,h);
%xr2 = x2(M+1:Nx-M);
xi2 = xh2(2*M+1:Nx);
%xa2 = xr2 + j*xi2;
%xh2 = abs(xa2);

t = (1/fs)*(0:1:length(x1i)-1);
%DSB-C
y_dsb_c1 = (1+x1i) .* cos(2*pi*fc1*t);
y_dsb_c2 = (1+x2i) .* cos(2*pi*fc2*t);
y_dsb_c = dA*(y_dsb_c1 + y_dsb_c2);


%DSB-SC
y_dsb_sc1 = (x1i) .* cos(2*pi*fc1*t);
y_dsb_sc2 = (x2i) .* cos(2*pi*fc2*t);
y_dsb_sc = dA*(y_dsb_sc1 + y_dsb_sc2);

% nowe t bo hilbert skróci³ 
t = (1/fs)*(0:1:length(xi1)-1);

%SSB-SC lewa '+'
x1i = x1i(M+1:Nx-M); %przycinam sygna³
x2i = x2i(M+1:Nx-M); %przycinam sygna³
y_ssb_sc_l1 = 0.5*x1i .* cos(2*pi*fc1*t) + 0.5*xi1 .* sin(2*pi*fc1*t);
y_ssb_sc_l2 = 0.5*x2i .* cos(2*pi*fc2*t) + 0.5*xi2 .* sin(2*pi*fc2*t);
y_ssb_sc_l = dA*(y_ssb_sc_l1 + y_ssb_sc_l2);

%SSB-SC prawa '-'
y_ssb_sc_p1 = 0.5*x1i .* cos(2*pi*fc1*t) - 0.5*xi1 .* sin(2*pi*fc1*t);
y_ssb_sc_p2 = 0.5*x2i .* cos(2*pi*fc2*t) - 0.5*xi2 .* sin(2*pi*fc2*t);
y_ssb_sc_p = dA*(y_ssb_sc_p1 + y_ssb_sc_p2);


%transformaty
HYdsb_c = fft(y_dsb_c)/length(y_dsb_c);
HYdsb_sc = fft(y_dsb_sc)/length(y_dsb_sc);
HYssb_sc1 = fft(y_ssb_sc_l)/length(y_ssb_sc_l);
HYssb_sc2 = fft(y_ssb_sc_p)/length(y_ssb_sc_p);

% Wykresy widm
f = (0:length(HYdsb_c)-1)/length(HYdsb_c)*fs;
figure;
subplot(1,2,1);
plot(f, abs(HYdsb_c));
title('fft DSB-C');
xlim([90e3 120e3]);
subplot(1,2,2);
plot(f, abs(HYdsb_sc));
title('fft DSB-SC');
xlim([90e3 120e3]);


f = (0:length(HYssb_sc1)-1)/length(HYssb_sc1)*fs;
figure;
subplot(211)
plot(f, abs(HYssb_sc1));
title('fft SSB-SC (+) lewa');
xlim([90e3 120e3]);
subplot(212);
plot(f, abs(HYssb_sc2));
title('fft SSB-SC (-) prawa');
xlim([90e3 120e3]);






