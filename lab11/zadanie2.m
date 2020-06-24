clc; clear all; close all;

[x,fpr] = audioread('DontWorryBeHappy.wav');
x = x(:,1);
N = 32; % d³ugoœæ okna
k = 0:N/2-1;
n = 0:N-1;
ileProbek = floor((length(x)-(N/2))/(N/2));
%okno
h = sin(pi*(n+0.5)/N);
A = sqrt(4/N) * cos((2*pi/N) * (k+0.5)' .* (n+0.5+N/4));
S = A.'; % transpozycja A
Q=2;
wynik=[];
for n = 1:ileProbek
    %analiza
    probka = x(1+(n-1)*(N/2) : N+(n-1)*(N/2));
    x1 = h' .* probka;
    x1 = A * x1;
    y = round(x1*Q);
    
    %synteza
    y = S * y;
    y = h .* y'; 
    y=y';
    wynik = [wynik;y]; % dla calego sygna³u
end

if(1)   % ³¹czenie sygna³u
  sygnal =[wynik(1:N/2)];
  for n=2:floor(length(wynik)/N)
      poprzednie = wynik(1+(n-2)*N : N+(n-2)*N);
      okno = wynik(1+(n-1)*N : N+(n-1)*N);
      suma = poprzednie(N/2+1:end) + okno(1:N/2);
      sygnal = [sygnal;suma];
  end
  sygnal = [sygnal;wynik(length(wynik)-N/2+1:end)];
 soundsc(sygnal,fpr);
end
