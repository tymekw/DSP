clc;
clear all;
close all;

fs = 256*10^3;% Hz
f0 = fs/2;
f = 0:1:fs;
w = 2*pi*f;
s = j*w;
N=6;
%% Butterworth
[b,a] = butter(N,2*pi*f0,'s');
%z = [j*2*pi*f0, -j*2*pi*f0];
p = roots(a);
z=roots(b);
%wzm2 = prod(-p);
%b = poly(z);

w = 0:0.1:20;

H = polyval(b,s) ./ polyval(a,s);
%H2 = H ./ max(H);

figure;
plot(real(z), imag(z), 'bo', real(p), imag(p), 'r*');
title("Zeros & Poles Butterworth"); grid;

figure;
subplot(211)
plot(f, 20*log10(abs(H)));
title("20*log10(|H(jw)|)"); xlabel("f[Hz]");
subplot(212);
plot(f, unwrap(angle(H))); grid;
title("faza [rad]"); xlabel("f[Hz]");

%% Czebyszew1

[b,a] = cheby1(N,3,2*pi*85*10^3,'s');
p = roots(a);
z=roots(b);

H = polyval(b,s) ./ polyval(a,s);

figure;
plot(real(z), imag(z), 'bo', real(p), imag(p), 'r*');
title("Zeros & Poles Czebyszew1"); grid;

figure;
subplot(211)
plot(f, 20*log10(abs(H)));
title("20*log10(|H(jw)|) Czebyszew1"); xlabel("f[Hz]");
subplot(212);
plot(f, unwrap(angle(H))); grid;
title("faza [rad]"); xlabel("f[Hz]");
%% Czebyszew2

[b,a] = cheby2(N,40,2*pi*f0,'s');
p = roots(a);
z=roots(b);

H = polyval(b,s) ./ polyval(a,s);

figure;
plot(real(z), imag(z), 'bo', real(p), imag(p), 'r*');
title("Zeros & Poles Czebyszew2"); grid;

figure;
subplot(211);
plot(f, 20*log10(abs(H)));
title("20*log10(|H(jw)|) Czebyszew2"); xlabel("f[Hz]");
subplot(212);
plot(f, unwrap(angle(H))); grid;
title("faza [rad]"); xlabel("f[Hz]");

%% eliptyczny

[b,a] = ellip(N,3,40,2*pi*120*10^3,'s');
p = roots(a);
z=roots(b);

H = polyval(b,s) ./ polyval(a,s);

figure;
plot(real(z), imag(z), 'bo', real(p), imag(p), 'r*');
title("Zeros & Poles eliptyczny"); grid;

figure;
subplot(211);
plot(f, 20*log10(abs(H)));
title("20*log10(|H(jw)|) eliptyczny"); xlabel("f[Hz]");
subplot(212);
plot(f, unwrap(angle(H))); grid;
title("faza [rad]"); xlabel("f[Hz]");

