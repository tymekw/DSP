clc;
clear all;

A = 230;  %amplitude   
T = 0.1;    %duration time
f0 = 50;    %signal frequency
fs3 = 200;    %sampling frequency
dt3 = 1/fs3;    %sampling period
t3 = 0:dt3:T;   %time moments for x(n)
x3 = A*sin(2*pi*f0*t3);    %generating sampled sin x(n)


dt = 1/10000;
t = 0:dt:T;
x = A*sin(2*pi*f0*t);    %generating analog sin x(n)

%x - analog signal(t) // x3 - sampled siganl (t3)
% t = nT =t3 where n = probe nbr and T = dt3

sum = zeros(size(t));
for n =1:length(t3)
    splot = x3(n)*sinc(fs3*(t-t3(n)));
    sum = sum+splot;
end

hold on
plot(t,x,'g-');%analogowy
plot(t3,x3,'b-');

plot(t,sum);%rekonstruowany
plot(t,abs(x-sum),'r-');%b³ad
legend("analog", "sampled","rekonstruowany","blad");
