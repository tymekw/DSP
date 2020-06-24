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
clear f;

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

%figure;
%psd(spectrum.welch('Hamming',1024), ym,'Fs',fs);

ym = ym(1:fs/bwSERV:end);
%soundsc(ym,bwSERV);

f=0:1:fpr/2;
H1 = freqz(h,1,f,fpr);

%{
figure;
hold all;
plot(f,20*log10(abs(H1)));
plot([0 60e3],[-3 -3],'r-'); % -3dB
plot([0 60e3],[-40 -40],'g-'); % -40dB
plot([15e3 15e3],[-300 50],'r-'); % 15kHz
plot([19e3 19e3],[-300 50],'g-'); % 19kHz
title('Filtr sk³adowej mono');
xlim([0 60e3]);
%}
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
fpilot = 19062;

f=0:1:fpr/2;
H1 = freqz(hp,1,f,fpr);

%{
figure;
hold all;
hold all;
plot(f,20*log10(abs(H1)),[0 60e3],[-3 -3],'r-');
plot([19e3 19e3],[-250 0],'r-'); % 19kHz
title('Filtr pilota');
xlim([0 60e3]);
%}



%% 2 podpunkt
f0 = 2*fpilot;
f1 = fpilot + 500; %Hz lewa strona
f2 = 3*fpilot - 500; %Hz prawa strona

MM = 1025; %d³ugoœæ filtru poprawna
M = (MM-1)/2;
n = -M:1:M; % filtr przechodz¹cy przez zero

%LP filter do fpilot Hz
hLP2p = 2*f2/fpr*sin(2*pi*f2/fpr*n)./(2*pi*f2/fpr*n);
hLP2p(M+1) = 2*f2/fpr;

%LP filter do 3* fpilot Hz
hLP1p = 2*f1/fpr*sin(2*pi*f1/fpr*n)./(2*pi*f1/fpr*n);
hLP1p(M+1) = 2*f1/fpr;

%BP filter
hb = hLP2p-hLP1p;

%BP + okna
hb = hb .* blackmanharris(MM)';


f=0:1:fpr/2;
H1 = freqz(hb,1,f,fpr);
figure;
hold all;
plot(f,20*log10(abs(H1)),'r-');
title('Filtr BP');
xlim([0 60e3]);

%% 3 podpunkt
%filtracja y = conv(x,h); y = filter(h,1,x);
clear wideband_signal;
clear wideband_signal_shifted;
clear wideband_signal_filtered;
clear hLP1;
clear hLP2;
clear hLP1p;
clear hLP2p;
clear dx;
clear H1;
clear a;
clear b;


Nx = length(y);
Nh = length(hb);
hx = zeros(1,Nh);
for k = 1:Nx
    hx = [y(k) hx(1:Nh-1)];
    ym1(k) = sum(hx .* hb);
end

clear y;



N = 1600000;
c = cos(2*pi*(fpilot/bwSERV)*[0:N-2]');
c=c';
ym1 = ym1 .* c;
figure;
psd(spectrum.welch('Hamming',1024), ym1,'Fs',fs);
title('Przesuniêcie do 0 Hz cos 2fpl');



%{
%% 4 podpunkt
%fpr = 30e3;




%{
L=500;
fs= fpr;
hLPaudio = fir1(500,(bwSERV/2)/(fs/2),kaiser(L+1,7));
ys = filter( hLPaudio, 1, ym );
% Pozostawienie co fs/Abw-tej próbki
  ys = ys(1:fs/bwSERV:end);
% Synchronizacja czasowa sygna³ów L+R i L-R (uwzglêdnienie opóŸnienia L-R)
  delay = (500/2)/(fs/bwSERV); ym = ym(1:end-delay); ys=2*ys(1+delay:end); % 2 od modul
  clear ymm yss;
% Odtworzenie kana³ów L i R
  y1 = 0.5*( ym + ys ); y2 = 0.5*( ym - ys ); clear ym ys;
% De-emfaza, ch-ka p³aska do 2.1 kHz, potem opadanie 20 dB na dekadê
  %y1 = filter(b_de,a_de,y1); y2 = filter(b_de,a_de,y2);
% Ods³uchanie sygna³u stereo po de-emfazie
  soundsc([y1' y2'],bwSERV); pause

%}


%% 4 podpunkt

N = 128;
b = fir1(N,[30/80e3 30e3/80e3],blackmanharris(N+1));
[h,f] = freqz(b,1,0.1:0.5:80e3,80e3);

y_stereo_filtered  = filtfilt(b,1,ym1);
%y_stereo_filtered = decimate(y_stereo_filtered,5);
y_stereo_filtered = resample(y_stereo_filtered,16,80);

%% widmo odfiltrowanego stereo do 30khz
figure(7)
psd(spectrum.welch('Hamming',1024), y_stereo_filtered,'Fs',fs/5);
title('Resampled Stereo');

%% yl and yr
ym = ym-mean(ym);
ym = ym/(1.001*max(abs(ym)));


y_stereo_filtered = y_stereo_filtered-mean(y_stereo_filtered);
y_stereo_filtered = y_stereo_filtered./(1.001*max(abs(y_stereo_filtered)));
yl = 0.5*(ym+y_stereo_filtered);
yr = 0.5*(ym-y_stereo_filtered);

y_kek = [yl yr];
%% Listen stereo to the final result 
soundsc( y_kek, bwAUDIO*2);
%% Listen mono to the final result 
%soundsc(ym,bwAUDIO*2)
%}