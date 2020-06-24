clear all; close all;

%wzorce kosinusowe DCT-II
%n -ta próbka k-tego sygna³u (funkcji bazowej)
%funcje (sygnaly) w wierszach

N=20;
s = sqrt(1/N); %s0
for k =1:20
    for n=1:20
        A(k,n) = s*cos(pi*(k-1)/N * ((n-1)+0.5));
    end
    s = sqrt(2/N); %sk
end

S = transpose(A);    %funkcjie (sygnaly) w kolumnach
I = S*A;             %czy jest to macierz jednostkowa???

I1 = eye(n);        %na pewno macierz jednostkowa

roznica = I-I1;
licznik = 1;
for i=1:N
    for j=1:N
        if( abs(roznica(i,j)) > (10^(-13)))
            fprintf("no niedobrze %u,%u \n",i,j);
        else
            fprintf("ok: %u\n",licznik);
            licznik = licznik + 1;
        end
    end
end






 
x = randn(20,1);   %sygna³ losowy (wektor)
X = A*x;           %analiza ( X -> ile jest kazdego)

rek = S*X;         %synteza

roznica2 =x-rek;   %blad miedzy x a rek 
licznik1=1;
for i=1:N
     if( abs(roznica2(i)) > (10^(-13)))
         fprintf("no niedobrze %u,%u \n",i,j);
     else
         fprintf("ok: %u\n",licznik1);         
     end
     licznik1 = licznik1 + 1;
end

%dla dociekliwych

