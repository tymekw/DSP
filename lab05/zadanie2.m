clc;
clear all;
close all;

f0=100;
w3db = 2*pi*f0; %rad/s
f= 0:1:200;
s= 1i * 2*pi*f;
%f =0:0.1:w3db/2*pi;

%% N=2
N=2;
for k = 1:N
        p2(k)= w3db*exp(1i*((pi/2)+(1/2)*(pi/N)+(k-1)*(pi/N)));
end
z2=[];
wzm2 = prod(-p2); %wzmocnienie to iloczyn zanegowanych biegunów

b2 = wzm2*poly(z2);
a2 = poly(p2);
H2 = polyval(b2,s) ./ polyval(a2,s);

%% N=4
N=4;
for k = 1:N
        p4(k)= w3db*exp(1i*((pi/2)+(1/2)*(pi/N)+(k-1)*(pi/N)));
end
z4=[];
wzm4 = prod(-p4); 

b4 = wzm4 * poly(z4);
a4 = poly(p4);
H4 = polyval(b4,s) ./ polyval(a4,s);

%% N=6
N=6;
for k = 1:N
        p6(k)= w3db*exp(1i*((pi/2)+(1/2)*(pi/N)+(k-1)*(pi/N)));
end
z6=[];
wzm6 = prod(-p6); 

b6 = wzm6 * poly(z6);
a6 = poly(p6);
H6 = polyval(b6,s) ./ polyval(a6,s);

%% N=8
N=8;
for k = 1:N
        p8(k)= w3db*exp(1i*((pi/2)+(1/2)*(pi/N)+(k-1)*(pi/N)));
end
z8=[];
wzm8 = prod(-p8); 

b8 = wzm8 * poly(z8);
a8 = poly(p8);
H8 = polyval(b8,s) ./ polyval(a8,s);

%% plots
figure;
hold on;
title("20*log10(|H(jw)|) liniowo");
plot(f,20*log10(abs(H2)));
plot(f,20*log10(abs(H4)));
plot(f,20*log10(abs(H6)));
plot(f,20*log10(abs(H8)));
legend("N=2","N=4","N=6","N=8");
grid;
xlabel("f [Hz]");
hold off;

figure;
title("20*log10(|H(jw)|) logarytmicznie");
semilogx(f,20*log10(abs(H2)));
hold on;
semilogx(f,20*log10(abs(H4)));
semilogx(f,20*log10(abs(H6)));
semilogx(f,20*log10(abs(H8)));
hold off;
legend("N=2","N=4","N=6","N=8");
grid;
xlabel("f [Hz]");

figure;
hold on;
title("angle H(jw) liniowo");
plot(f,unwrap(angle(H2)));
plot(f,unwrap(angle(H4)));
plot(f,unwrap(angle(H6)));
plot(f,unwrap(angle(H8)));
legend("N=2","N=4","N=6","N=8");
grid;
xlabel("f [Hz]");
hold off;

%% odpowiedz impulsowa i skok jednostkowy
%H = tf(a4,b4);
%printsys(a4,b4,'s') %'s' bo analog

figure;
impulse(b4,a4);
figure;
step(b4,a4);

