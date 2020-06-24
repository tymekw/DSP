clc;
clear all;
close all;

fs= 3.2e6; %sampling freq
N = 32e6; %number of samples
fc = 0.40e6; %central freq MF station

bwSERV = 80e3; %  % bandwidth of an FM service (bandwidth ~= sampling frequency!)
bwAUDIO = 16e3;     % bandwidth of an FM audio (bandwidth == 1/2 * sampling frequency!)

f = fopen('samples_100MHz_fs3200kHz.raw');
s = fread(f, 2*N, 'uint8');
fclose(f);

s=s-127;
%% IQ -> complex
wideband_signal = s(1:2:end) + sqrt(-1)*s(2:2:end); clear s;

%% Extract carrier of selected service, then shift in frequency the selected service to the baseband
wideband_signal_shifted = wideband_signal .* exp(-1i*2*pi*fc/fs*[0:N-1]');

%% Filter out the service from the wide-band signal
%Lab6
[b,a] = butter(4,bwSERV/fs);
wideband_signal_filtered = filter( b, a, wideband_signal_shifted );

%% Down-sample to service bandwidth - bwSERV = new sampling rate
x = decimate(wideband_signal_filtered,20);
fs=fs/20;

%% FM demodulation
dx = x(2:end).*conj(x(1:end-1));
y = atan2( imag(dx), real(dx) );

figure;
psd(spectrum.welch('Hamming',1024), y,'Fs',fs);

%% filtering mono

fpr = fs; %Hz
f1 = 30; %Hz lewa strona
f2 = 15e3; %Hz prawa strona
%MM = 128; %d³ugoœæ flitru z³a bo parzysta
MM = 129; %d³ugoœæ filtru poprawna
M = (MM-1)/2;
n = -M:1:M; % filtr przechodz¹cy przez zero

%LP filter do 15e3Hz
hLP2 = 2*f2/fpr*sin(2*pi*f2/fpr*n)./(2*pi*f2/fpr*n);
hLP2(M+1) = 2*f2/fpr;

%LP filter do 30Hz
hLP1 = 2*f1/fpr*sin(2*pi*f1/fpr*n)./(2*pi*f1/fpr*n);
hLP1(M+1) = 2*f1/fpr;

%BP filter
h = hLP2-hLP1;

%BP + okna
h = h .* chebwin(129)';


Nx = length(y);
Nh = length(h);
hx = zeros(1,Nh);
for k = 1:Nx
    hx = [y(k) hx(1:Nh-1)];
    ym(k) = sum(hx .* h);
end

figure;
psd(spectrum.welch('Hamming',1024), ym,'Fs',fs);

ym = ym(1:fs/bwSERV:end);
%soundsc(ym,bwSERV);

f=0:1:fpr/2;
H1 = freqz(h,1,f,fpr);

figure;
hold all;
plot(f,20*log10(abs(H1)));
plot([0 60e3],[-3 -3],'r-'); % -3dB
plot([0 60e3],[-40 -40],'g-'); % -40dB
plot([15e3 15e3],[-300 50],'r-'); % 15kHz
plot([19e3 19e3],[-300 50],'g-'); % 19kHz
title('Filtr sk³adowej mono');
xlim([0 60e3]);

%% pilot filtering

fpr = fs; %Hz
f0 = 19e3;
f1 = f0-400; %Hz lewa strona
f2 = f0+400; %Hz prawa strona
%MM = 128; %d³ugoœæ flitru z³a bo parzysta
MM = 501; %d³ugoœæ filtru poprawna
M = (MM-1)/2;
n = -M:1:M; % filtr przechodz¹cy przez zero

%LP filter do 19300Hz
hLP2p = 2*f2/fpr*sin(2*pi*f2/fpr*n)./(2*pi*f2/fpr*n);
hLP2p(M+1) = 2*f2/fpr;

%LP filter do 18700Hz
hLP1p = 2*f1/fpr*sin(2*pi*f1/fpr*n)./(2*pi*f1/fpr*n);
hLP1p(M+1) = 2*f1/fpr;

%BP filter
hp = hLP2p-hLP1p;

%BP + okna
hp = hp .* chebwin(MM)';

Nx = length(y);
Nh = length(hp);
hx = zeros(1,Nh);
for k = 1:Nx
    hx = [y(k) hx(1:Nh-1)];
    ymp(k) = sum(hx .* hp);
end

figure;
psd(spectrum.welch('Hamming',1024), ymp,'Fs',fs);

%ym = ym(1:fs/bwSERV:end);
%soundsc(ym,bwSERV);

f=0:1:fpr/2;
H1 = freqz(hp,1,f,fpr);

figure;
hold all;
hold all;
plot(f,20*log10(abs(H1)),[0 60e3],[-3 -3],'r-');
plot([19e3 19e3],[-250 0],'r-'); % 19kHz
title('Filtr pilota');
xlim([0 60e3]);

figure;
spectrogram( ymp, 4096, 4096-512, [0:5:20000], 2*bwSERV );

%{
  fcentr = 19e3; df1 = 1000; df2 = 2000;
  ff = [ 0 fcentr-df2 fcentr-df1 fcentr+df1 fcentr+df2 fs/2 ]/(fs/2);
  fa = [ 0 0.01 1 1 0.01 0 ];
  hBP19 = firpm(500,ff,fa);

p = filter(hBP19,1,y);


figure;
psd(spectrum.welch('Hamming',1024), p,'Fs',fs);

%ym = ym(1:fs/bwSERV:end);
%soundsc(ym,bwSERV);

f=0:1:fpr/2;
H1 = freqz( hBP19 ,1,f,fpr);

figure;
hold all;
hold all;
plot(f,20*log10(abs(H1)),[0 60e3],[-3 -3],'r-');
plot([19e3 19e3],[-250 0],'r-'); % 19kHz
title('Filtr pilota');
xlim([0 60e3]);
%}











