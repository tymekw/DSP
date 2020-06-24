clc; clear all; close all;
[x,fpr] = audioread('DontWorryBeHappy.wav', 'native');
x = double(x);
%soundsc(x,fpr);
a = 0.9545;
Nb = 8; %liczba bitów %od 8 brzmi ok
d = x - a*[[0,0]; x(1:end-1,:)];
dq = lab11_kwant(d,Nb);

% dekoder
% koder
% d(n) = x(n) – a*x(n-1);
% y(n) =dq(n) +a*y(n-1);

%dekoder
yd=[];
for k=1:2
   dq1 = dq(:,k);
   y(1) = dq1(1);
   for n=2:length(dq)
      y(n) = dq1(n) + a*y(n-1);  
   end
   yd(:,k) = y;
end

% porównanie sygna³ów
figure;
hold on;
plot(x(:,1),'b');
plot(yd(:,1),'r')
legend("oryginalny","zrekonstruowany");
title('lewy kana³')

figure;
hold on;
plot(x(:,2),'b');
plot(yd(:,2),'r')
legend("oryginalny","zrekonstruowany");
title('prawy kana³')
%soundsc(yd,fpr);


figure;
plot(abs(x(:,2).^2-yd(:,2).^2));


figure;
n=1:length(x);
plot(n,x,'b',n,dq,'r');