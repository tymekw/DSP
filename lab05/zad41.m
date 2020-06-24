clear all; close all; clc;
 
N = 4;
 
freq_req = 96000000;
 
% limits : 95900000 && 96100000 +/- 100 kHz
% limits : 95000000 && 97000000 +/- 1 MHz
 
fs = freq_req*2; % Hz
Ws = 2*pi * fs;
 
f = 0:fs+2e6;
w = 2*pi * f;
s = 1i * 2*pi * f;
 


% [b_band, a_band] = cheby1(N, 3, [Ws/2-1e5 Ws/2+1e5], 's');
% [b_band, a_band] = cheby2(N, 40, [Ws/2-1e6 Ws/2+1e6], 's');
[b_band, a_band] = ellip(N, 3, 40, [Ws/2-6e5 Ws/2+6e5], 's');
 
f_band = f( freq_req-2e6 : freq_req+2e6 );
H_band = freqs(b_band, a_band, 2*pi*f_band);
 
figure;
plot(f_band, 20*log10( abs( H_band ) ), 'b-' );
title("Bandpass filter RMF");
xlabel("Frequency [ Hz ]"); ylabel("Magnitude [ dB ]");
xlim([96e6-2e6 96e6+2e6]); ylim([-60 20]);
yline(-40, 'r--'); yline(-3, 'r--'); yline(0, 'r--');
xline(95e6, '--m'); xline(97e6, '--m'); % +/- 1 MHz
xline(961e5, '--k'); xline(959e5, '--k'); % +/- 100 kHz
grid;