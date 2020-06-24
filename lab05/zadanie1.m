clc;
clear all;
close all;

p1 = -0.5 + 1i*9.5;
p2 = -0.5 - 1i*9.5;
p3 = -1 + 1i*10;
p4 = -1 - 1i*10;
p5 = -0.5 + 1i*10.5;
p6 = -0.5 - 1i*10.5;
z1 = 1i*5;
z2 = -1i*5;
z3 = 1i*15;
z4 = -1i*15;

p = [p1,p2,p3,p4,p5,p6];    % bieguny
z = [z1,z2,z3,z4];          % zera s¹ na osi jw wiêc sie wyzeruje
b = poly(z);
a = poly(p);


f = 0:0.01:3;
w = 2*pi*f;
%w = 0:0.1:20;
s = j*w;
H = polyval(b,s) ./ polyval(a,s);

figure;
plot(real(z), imag(z), 'bo', real(p), imag(p), 'r*');
legend("zeros", "poles"); title("Zeros & Poles"); grid;

figure;
plot(w, 20*log10(abs(H)));
title("20*log10(|H(jw)|)"); xlabel("w[rad/s]");

figure;
plot(w, abs(H));
title("|H(jw)|"); xlabel("w[rad/s]");

%modyfikacja wzmocnienia
H2 = H ./ max(H);

figure;
plot(w, abs(H2));
title("|H(jw)| wzmocnienie 1"); xlabel("w[rad/s]");

figure;
plot(w, 20*log10(abs(H2)));
title("20*log10(|H(jw)|) wzmocnienie 1"); xlabel("w[rad/s]");


%ch-f-cz
figure;
plot(w, unwrap(angle(H))); grid;
title("faza [rad]"); xlabel("w[rad/s]");

