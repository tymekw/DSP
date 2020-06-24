clc;
clear all;

% czas trwania 10s
fs = 44100; %czestotliwosc probkowania
dt = 1/fs;
tStart = 0;
tEnd = 10;
t = tStart:dt:tEnd;
N1  = ceil(10/dt);

f1 = 1000;  %czestotliwosc poczatkowa
fd = 5000;  %zmiana /1s

y = chirp(t,f1,1,f1+fd,'linear',90);
plot(t,y);
sound(y,fs);