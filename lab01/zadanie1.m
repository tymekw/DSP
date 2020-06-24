clc;
clear all;

A=230;
f=50;
T = 0.1;

fs1=10000;
fs2=500;
fs3=200;

dt1=1/fs1;
dt2=1/fs2;
dt3=1/fs3;

t1=0:dt1:T;
t2=0:dt2:T;
t3=0:dt3:T;

x1 = A*sin(2*pi*f*t1);
x2 = A*sin(2*pi*f*t2);
x3 = A*sin(2*pi*f*t3);

figure;
plot(t1,x1,'b-');
hold on;
plot(t2,x2,'r-o');
plot(t3,x3,'k-x');
legend("10kHz","500Hz","200Hz");

%---------------------------------------

T=1;
f=50;

fs1=10000;
fs2=26;
fs3=25;
fs4=24;

dt1=1/fs1;
dt2=1/fs2;
dt3=1/fs3;
dt4=1/fs4;

t1=0:dt1:T;
t2=0:dt2:T;
t3=0:dt3:T;
t4=0:dt4:T;

x1 = A*sin(2*pi*f*t1);
x2 = A*sin(2*pi*f*t2);
x3 = A*sin(2*pi*f*t3);
x4 = A*sin(2*pi*f*t4);

figure;
plot(t1,x1,'b-');
hold on;
plot(t2,x2,'g-o');
plot(t3,x3,'r-o');
plot(t4,x4,'k-o');
xlim([0 1]);
legend("10kHz","26kHz","25kHz","24kHz");

%c-----------------------------------------------
clc;
clear all;

T=1;
fs=100;
dt=1/fs;
t=0:dt:T;

f=0;
x=zeros([61 101]);
figure;
for i=1:61
    y = sin(2*pi*t*f);
    for j=1:101
     x(i,j)=y(j);  
    end
    disp("numer obiegu: ");
    disp(i);
    disp("wartoœæ czêstotliwoœci: ");
    disp(f);
    plot(t,y);
    pause;
    f = f+5;
end

figure;
subplot(3,1,1);
hold all;
plot(t,x(2,:), 'g-o');
plot(t,x(22,:), 'r-x');
plot(t,x(42,:), 'k-o');
title('Próbkowanie sinusa 50Hz fs= 5 105 205');
legend('5Hz','105Hz','205Hz');
xlabel('Czas [s]');

subplot(3,1,2);
hold all;
plot(t,x(20,:), 'g-o');
plot(t,x(40,:), 'r-x');
plot(t,x(60,:), 'k-o');
title('Próbkowanie sinusa 50Hz Fs = 95 195 295');
legend('95Hz','195Hz','295Hz');
xlabel('Czas [s]');

subplot(3,1,3);
hold all;
plot(t,x(20,:), 'g-o');
plot(t,x(22,:), 'r-o');
xlim([0 1]); %skalowanie osi x - czas%
title('Próbkowanie sinusa 50Hz Fs = 95 105');
legend('95Hz','105Hz');
xlabel('Czas [s]');

    
%dla cos
clc;
clear all;

T=1;
fs=100;
dt=1/fs;
t=0:dt:T;

f=0;
x=zeros([61 101]);
figure;
for i=1:61
    y = cos(2*pi*t*f);
    for j=1:101
     x(i,j)=y(j);  
    end
    disp("numer obiegu: ");
    disp(i);
    disp("wartoœæ czêstotliwoœci: ");
    disp(f);
    plot(t,y);
    pause;
    f = f+5;
end

figure;
subplot(3,1,1);
hold all;
plot(t,x(2,:), 'g-o');
plot(t,x(22,:), 'r-x');
plot(t,x(42,:), 'k-o');
%xlim([0 1]); %skalowanie osi x - czas%
title('Próbkowanie cosinusa 50Hz fs= 5 105 205');
legend('5Hz','105Hz','205Hz');
xlabel('Czas [s]');

subplot(3,1,2);
hold all;
plot(t,x(20,:), 'g-o');
plot(t,x(40,:), 'r-x');
plot(t,x(60,:), 'k-o');
%xlim([0 1]); %skalowanie osi x - czas%
title('Próbkowanie cosinusa 50Hz Fs = 95 195 295');
legend('95Hz','195Hz','295Hz');
xlabel('Czas [s]');


subplot(3,1,3);
hold all;
plot(t,x(20,:), 'g-o');
plot(t,x(22,:), 'r-o');
xlim([0 1]); %skalowanie osi x - czas%
title('Próbkowanie cosinusa 50Hz Fs = 95 105');
legend('95Hz','105Hz');
xlabel('Czas [s]');












