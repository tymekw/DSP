clear all;
clc;

N=100;          %liczba probek
fs = 1000;      %czestotliwosc probkowania Hz
t=(0:N-1)*1/fs;

%czestotliwosci sinusoid
f1=50;
f2=100;
f3=150;

%amplitudy sinusoid;
A1=50;
A2=100;
A3=150;

%sygnaly skladowe
x1 = A1*cos(2*pi*f1*t);
x2 = A2*cos(2*pi*f2*t);
x3 = A3*cos(2*pi*f3*t);

x = x1+x2+x3;


%macierz A=DCT
%n -ta próbka k-tego sygna³u (funkcji bazowej)
%funcje (sygnaly) w wierszach

s = sqrt(1/N);
for k=1:N
    for n=1:N
        A(k,n)=s*cos(pi*(k-1)/N *((n-1)+0.5));
    end
    s = sqrt(2/N);
end


%maciesz S=IDCT
S = A';

%wyœwietlanie wartosci wierszy A i kolumn S
%{
figure;
for i=1:N
    fprintf("\nWIERSZ %u DCT  |  KOLUMNA %u IDCT \n\n",i,i);
    for j =1:N
        fprintf("%u  |  %u \n",A(i,j),S(j,i));
    end
    plot(A(i,:));
    hold on; 
    plot(S(:,i));
    legend("A","S");
    pause;
    hold off;
end
%}

%analiza y=Ax

y = A*x';   %y => ile jest danych czest z A w x

disp("Wartoœci y(1:N)");
disp(y(1:N));

%wyci¹gniêcie wspó³czynników niezerowych
k=1;
numery = zeros(1,100);
for n=1:100
    if(abs(y(n))>(10^(-11)))
        wspol_niezer(k) = y(n);
        numery(k) = n;
        k = k+1;
    end
end

figure;
subplot(221);
stem(wspol_niezer);
title("wart. wspó³.");
subplot(223);
stem(x);
title("wart. amplitud");

f = (0:N-1)*500/N;
subplot(222);
stem(f,y*sqrt(2/N));
title("wart. wspó³.");
subplot(224);
stem(f,x);
title("wart. amplitud");


%sprawdzenie rekonstrukcji

xr = S*y;   %rekonstrukcja

roznica =x-xr';   %blad miedzy x a xr 
licznik1=1;
for i=1:N
     if( abs(roznica(i)) > (10^(-11)))
         fprintf("no niedobrze %u \n",i);
     else
         fprintf("ok: %u \n",licznik1);         
     end
     licznik1 = licznik1 + 1;
end

maxblad = max(abs(roznica));
fprintf("max blad wynosi: %d \n", maxblad);

%x_1 zmiana czestotliwosci f2 na 105
f2=105;
x2 = A2*cos(2*pi*f2*t);
x_1=x1+x2+x3;

f = (0:N-1)*500/N;


y_1=A*x_1';
figure;
hold on;
stem(f,y_1*sqrt(2/N),'b');
stem(f,y*sqrt(2/N),'r');
legend("przesuniete", "normalne");
title("f2 przesuniete do 105Hz");
%sprawdzenie rekonstrukcji

xr1 = S*y_1;   %rekonstrukcja

roznica1 =x_1-xr1';   %blad miedzy x a xr 
licznik1=1;
for i=1:N
     if( abs(roznica1(i)) > (10^(-11)))
         %fprintf("no niedobrze %u \n",i);
     else
         %fprintf("ok: %u \n",licznik1);         
     end
     licznik1 = licznik1 + 1;
end

maxblad1 = max(abs(roznica1));
fprintf("max blad wynosi: %d \n", maxblad1);






%x_2 zmiana czestotliwosci o 2.5 
f1=52.5;
f2=102.5;
f3=152.5;
x1 = A1*cos(2*pi*f1*t);
x3 = A3*cos(2*pi*f3*t);
x2 = A2*cos(2*pi*f2*t);

x_2=x1+x2+x3;

y_2=A*x_2';
figure;
stem(f,y_2*sqrt(2/N),'b');
hold on;
stem(f,y*sqrt(2/N),'r');
legend("przesuniete", "normalne");
title("wszytsko przesuniete o 2.5Hz");
%sprawdzenie rekonstrukcji

xr2 = S*y_2;   %rekonstrukcja

roznica2 =x_2-xr2';   %blad miedzy x a xr 
licznik2=1;
for i=1:N
     if( abs(roznica2(i)) > (10^(-11)))
         %fprintf("no niedobrze %u \n",i);
     else
         %fprintf("ok: %u \n",licznik1);         
     end
     licznik1 = licznik1 + 1;
end

maxblad2 = max(abs(roznica2));
fprintf("max blad wynosi: %d \n", maxblad2);

figure;
plot(f,x,'b');
hold on;
plot(f,xr1,'r');
plot(f,xr2,'g');
