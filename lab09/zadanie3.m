clc;
clear all;
close all;

%udawany pilot 1
fpilot = 19e3; %czest pilota
fs = 40000; %czest probk
dt = 1/fs;
T=1;
t = 0:dt:T-dt; %czas trwania sygna³u
con_phase = pi/4;
p = cos(2*pi*fpilot*t + con_phase);


% Petla PLL z filtrem typu IIR do odtworzenia czêstotliwoœci i fazy pilota [7]
% i na tej podstawie sygna³ów noœnych: symboli c1, stereo c38 i danych RDS c57
freq = 2*pi*fpilot/fs;
theta = zeros(1,length(p)+1);
alpha = 1e-2;
beta = alpha^2/4;
for n = 1 : length(p)
perr = -p(n)*sin(theta(n));
theta(n+1) = theta(n) + freq + alpha*perr;
freq = freq + beta*perr;
end

blad = 2*pi*fpilot*t+pi/4 - theta(1:end-1);

figure
plot(t, blad); title('1 podpunkt');



%udawany pilot 2
df = 10;
fm = 0.1;
fpilot = 19e3; %czest pilota
fs = 8e3; %czest probk
dt = 1/fs;
T=1;
t = 0:dt:T-dt; %czas trwania sygna³u
phase = df*sin(1*pi*fm*t);
p = cos(2*pi*fpilot*t + phase);


% Petla PLL z filtrem typu IIR do odtworzenia czêstotliwoœci i fazy pilota [7]
% i na tej podstawie sygna³ów noœnych: symboli c1, stereo c38 i danych RDS c57
freq = 2*pi*fpilot/fs;
theta = zeros(1,length(p)+1);
alpha = 1e-2;
beta = alpha^2/4;
for n = 1 : length(p)
perr = -p(n)*sin(theta(n));
theta(n+1) = theta(n) + freq + alpha*perr;
freq = freq + beta*perr;
end

blad = 2*pi*fpilot*t+phase- theta(1:end-1);
figure
plot(t, blad); title('2 podpunkt');




%udawany pilot 3
fpilot = 19e3; %czest pilota
fs = 8e3; %czest probk
dt = 1/fs;
T=1;
t = 0:dt:T-dt; %czas trwania sygna³u
phase = pi/4;
p1 = cos(2*pi*fpilot*t + phase);
%p = awgn(p1,0,'measured');
%p = awgn(p1,5,'measured');
%p = awgn(p1,10,'measured');
p = awgn(p1,20,'measured');


% Petla PLL z filtrem typu IIR do odtworzenia czêstotliwoœci i fazy pilota [7]
% i na tej podstawie sygna³ów noœnych: symboli c1, stereo c38 i danych RDS c57
freq = 2*pi*fpilot/fs;
theta = zeros(1,length(p)+1);
alpha = 1e-2;
beta = alpha^2/4;
for n = 1 : length(p)
perr = -p(n)*sin(theta(n));
theta(n+1) = theta(n) + freq + alpha*perr;
freq = freq + beta*perr;
end

blad = abs(p- cos(theta(1:end-1)));
figure
plot(t, blad); title('3 podpunkt');

for i=1:length(blad)
    if(abs(blad(i)) < 1e-4)
        display(i)
        break
    end
end

