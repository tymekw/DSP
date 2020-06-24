clear all; close all;

fs = 3.2e6;         % sampling frequency
N  = 32e6;         % number of samples (IQ)
fc = 0.40e6;        % central frequency of MF station

bwSERV = 80e3;     % bandwidth of an FM service (bandwidth ~= sampling frequency!)
bwAUDIO = 16e3;     % bandwidth of an FM audio (bandwidth == 1/2 * sampling frequency!)

f = fopen('samples_100MHz_fs3200kHz.raw');
s = fread(f, 2*N, 'uint8');
fclose(f);

s = s-127;

% IQ --> complex
wideband_signal = s(1:2:end) + sqrt(-1)*s(2:2:end); clear s;

% Extract carrier of selected service, then shift in frequency the selected service to the baseband
wideband_signal_shifted = wideband_signal .* exp(-1i*2*pi*fc/fs*[0:N-1]');

% Filter out the service from the wide-band signal

%%%%%%%%%

% A_stop = 40;		% t³umienie w pierwszym pa�mie zaoprowym
% F_stop = 80e3+2000;		% granica pasma zaporowego
% F_pass = 80e3;	% granica pasma przepustowego
% A_pass = 3;		% zafalowania w pa�mie przepustowym nie mog¹ byæ wiêksze ni¿ 3 dB
% 
% LowPassSpecObj1 = fdesign.lowpass(F_pass, F_stop, A_pass, A_stop,3.2e6);
% LowPassFilt1 = design(LowPassSpecObj1, 'ellip');
% 
% 
% [b,a] = tf(LowPassFilt1);
[b,a] = butter(4,80e3/1.6e6);

%%%%%%%%%

wideband_signal_filtered = filter( b, a, wideband_signal_shifted );
%wideband_signal_filtered = filter( LowPassFilt1, wideband_signal_shifted );

% Down-sample to service bandwidth - bwSERV = new sampling rate
%x = wideband_signal_filtered( 1:fs/(bwSERV*2):end );
x = decimate(wideband_signal_filtered,20);
% FM demodulation
dx = x(2:end).*conj(x(1:end-1));
y = atan2( imag(dx), real(dx) );

% Decimate to audio signal bandwidth bwAUDIO

%%%%%%%%%

f_norm=bwAUDIO/bwSERV;
[ba,aa]=butter(4,f_norm);
y = filter( ba, aa, y ); % antyaliasing filter
ym = decimate(y,5);
%ym = y(bwSERV/bwAUDIO);  % decimate (1/5)

%%%%%%%%%

% De-emfaza
% (...)

%f_norm=2.1e3/bwAUDIO;
%[baa,aaa]=butter(1,f_norm);
%ym = filter( baa, aaa, ym );

f_norm=2.1e3/bwAUDIO;
[bd,ad]=butter(1,f_norm);
ym = filter(bd,ad,ym);

%ym = filter(1,[1 -0.9735],ym);

% Listen to the final result
ym = ym-mean(ym);
ym = ym./(1.001*max(abs(ym)));
soundsc( ym, bwAUDIO*2);

psd(spectrum.welch('Hamming',1024), wideband_signal(1:N),'Fs',fs);
plot(abs(wideband_signal));

