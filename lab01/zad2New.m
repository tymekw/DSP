clc;
clear all;

% Program parameters
f0=50;  %freq of x(n)
fs = 200;   %freq of sampling
dt = 1/fs; %okres probkowania
czas = 0.1; %czas trwania sygnalu
A = 230;    %amplituda

%generation of sine x(n)  with f0 freq


t = 0:dt:czas;  %momenty gdy pobieramy probki
x = A * sin(2*pi*f0*t); %probki
N = length(x);      %liczba probek
figure
plot(t,x,'b-', t,x,'ro');
grid;
title("Sinus próbkowany z wysok¹ czêstotliwoœci¹");

%sinc(a) sin(a)/a function
y=zeros(1,N);
%for k=1:length(xO);

%odtwarzanie
fs1 = 10000; %czestotliwosc jak chce miec dokladnie odtworzone
T = 1/fs1;   %okres probkowania
ts = 0:T:czas;  %momenty gdy pobieramy probki
NX = length(ts); %ilosc probek szukanych

xAnalog = A*sin(2*pi*f0*ts);



for k=1:NX
    t1 = ts(k);
    y(k)=0;
    for n=1:N
        y(k) = y(k) + x(n)*sinc(fs*(ts(k)-t(n)));
    end
end
figure;
hold on;
plot(t,x,'b-');

plot(ts,y);
plot(ts,abs(y-xAnalog));
plot(ts,xAnalog);
legend("sprobkowany","zrekonstruowany","blad","analogowy");  
    











