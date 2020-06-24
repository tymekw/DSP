clear all; close all; clc;
[x fs] = audioread("mowa.wav");%jak jest tylko jedna kolumna to jest tak samo na prawym i lewym
plot(x);

n = [1112 4002 6462 8041 13740 15470 17030 19300 26400 35280];
A = dctmtx(256);
figure(2);
for i=1:length(n)
    Xk(i,1:256) = x(n(i):n(i)+255);
    b = x(n(i):n(i)+255);
    Yk(i,1:256) = A*x(n(i):n(i)+255);
    
    subplot(2,1,1);
    plot(Xk(i,:)); title("wartosci probek");
    subplot(2,1,2);
    plot((1:256)*fs/256, Yk(i,:)); title("wartosci wspolczynnikow");
    pause();
end
z = zeros(1,length(x));
y = [x z'];
% sound(y,fs);