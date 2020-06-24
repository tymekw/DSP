close all;
clear all;

fs = 3.2e6;         % sampling frequency
N  = 32e6;         % number of samples (IQ)
fc = 0.40e6;        % central frequency of MF station

bwSERV = 80e3;     % bandwidth of an FM service (bandwidth ~= sampling frequency!)
bwAUDIO = 16e3;     % bandwidth of an FM audio (bandwidth == 1/2 * sampling frequency!)

f = fopen('samples_100MHz_fs3200kHz.raw');
s = fread(f, 2*N, 'uint8');
fclose(f);

s = s-127;

%% IQ --> complex
wideband_signal = s(1:2:end) + sqrt(-1)*s(2:2:end); clear s;

%% Extract carrier of selected service, then shift in frequency the selected service to the baseband
wideband_signal_shifted = wideband_signal .* exp(-sqrt(-1)*2*pi*fc/fs*[0:N-1]');

%% Filter out the service from the wide-band signal
%Lab6
[b,a] = butter(4,bwSERV/fs);
wideband_signal_filtered = filter( b, a, wideband_signal_shifted );

%% Down-sample to service bandwidth - bwSERV = new sampling rate
x = decimate(wideband_signal_filtered,20);
fs = fs/20;

%% FM demodulation
dx = x(2:end).*conj(x(1:end-1));
y = atan2( imag(dx), real(dx) );

figure;
psd(spectrum.welch('Hamming',1024), y,'Fs',fs);
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% filtr sk³adowej mono 30Hz - 15kHz
N=128;
b = fir1(N,[30/bwSERV 15e3/bwSERV],blackmanharris(N+1));
[h,f] = freqz(b,1,0.1:0.05:80e3,80e3);


figure(1);
subplot(2,1,1);
hold all;
plot(2*f,20*log10(abs(h)));
plot([0 60e3],[-3 -3],'r-'); % -3dB
plot([0 60e3],[-40 -40],'g-'); % -40dB
plot([15e3 15e3],[-300 50],'r-'); % 15kHz
plot([19e3 19e3],[-300 50],'g-'); % 19kHz
title('Filtr sk³adowej mono');
xlim([0 60e3]);

y_audio =  filtfilt(b,1,y);
subplot(2,1,2);
psd(spectrum.welch('Hamming',1024), y_audio,'Fs',fs);
title('Przefiltrowane od 30Hz do 15kHz')

%% filtr pilota 19kHz +/-10Hz
N=512;
b = fir1(N,[(18.990e3)/bwSERV (19.010e3)/bwSERV],blackmanharris(N+1));
[h,f] = freqz(b,1,0.1:0.05:80e3,80e3);

figure(2);
subplot(2,1,1);
hold all;
plot(2*f,20*log10(abs(h)),[0 60e3],[-3 -3],'r-');
plot([19e3 19e3],[-250 0],'r-'); % 19kHz
title('Filtr pilota');
xlim([0 60e3]);

y_pik =  filtfilt(b,1,y);
%19,062 kHz
subplot(2,1,2);
psd(spectrum.welch('Hamming',1024), y_pik,'Fs',fs);
title('Pilot');


figure;
spectrogram( y_pik, 4096, 4096-512, [0:5:20000], 2*bwSERV );
