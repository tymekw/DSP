clc; clear all; close all;

x = open("lab08_am.mat");
x = x.s2;
Nx = length(x);
T=1; fs=1000; fc=200;
N=T*fs; dt = 1/fs;
t = (0:fs*T -1)*dt;
%filtr hilberta
M=60; n =-M:M; MM = 2*M+1;
hH = (1-cos(pi*n))./(pi*n); hH(M+1)=0;
h=hH;


%odpowiedŸ impulsowa filtru
f=-fs/2:1:fs/2;
H = freqz(h,1,f,fs);


%figure;
%stem(h); title("h(n) wagi filtra"); grid;

%w = chebwin(MM,90)';
%w=blackman(2*M+1)';
w=kaiser(2*M+1,10)';
h=h.*w;




%odpowiedŸ impulsowa filtru z oknem
H1 = exp(j*2*pi*f/fs*M).* freqz(h,1,f,fs); %exp zeby byl k¹t f/fs omega M przesuniecie
%figure;
%stem(abs(h)); title("hw(n) wagi filtra z oknem"); grid;
figure;
hold on;
plot(f,20*log10(abs(H)));
plot(f,20*log10(abs(H1)));
xlabel('f [Hz]'); title("|H(f)|"); grid;
legend("zwyk³y","okno");

%{
figure;
plot(f,pi-unwrap(angle(H1)));
title("K¹t");
%}


% filtracja
y = conv(x,h);
%{
Nh = length(h);
bx = zeros(1,Nh);
for n=1:N
    bx = [x(n) bx(1:Nh-1)];
    y(n) = sum(bx .* h);
end
%}

%sk³adam sygna³ y z x -> tworzê analityczny
xr = x(M+1:Nx-M);
xi = y(2*M+1:Nx);

xa = xr + j*xi;
AM = abs(xa);

%AM == m
m = sqrt(xr.^2 + xi.^2);
%m =AM;
figure;
hold on;
plot(AM); plot(xr); plot(xi);
legend("sygna³ moduluj¹cy m","sygna³ wejœciowy x","sygna³ przesuniêty Hilbertem xh");

%analiza czêstotliwoœciowa
mFFT = abs(fft(m))/length(m);
mFFT(1) = 0;
f=fs*(0:length(m)-1)/length(m);
figure;
plot(f,mFFT);

A1=0.36;%0.24;
A2=0.07;%0.1;
A3=0.25;%0.15;
f1=9;%2;
f2=30;%40;
f3=60;%60;


%t=t(M+1:end-M);
mt=1+2*A1*cos(2*pi*f1*t)+2*A2*cos(2*pi*f2*t)+2*A3*cos(2*pi*f3*t);
mt = mt(M+1:end-M);
%mt=abs(mt);
figure;
plot(m);
hold on;
plot(mt);
legend("obliczony" , "odtworzony");
title("obliczony i odtworzony");


mFFT = abs(fft(mt))/length(mt);
mFFT(1) =0;
f=fs*(0:length(mt)-1)/length(mt);
figure;
plot(f,mFFT);
