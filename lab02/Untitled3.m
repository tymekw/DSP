clear all;
clc;

N = 100; %liczba probek%
fs = 1000; %czestotliwosc probkowania%
T = 0.1; %czas trwania probkowania (100 probek dla 1000Hz = 0.1s)%

% Czestotliwosci sinusoid %
f1 = 50; 
f2 = 100;
f3 = 150;

% Okresy pr?bkowa? %
fs1 = 1/fs;  
fs2 = 1/fs;    
fs3 = 1/fs;  

% Pr?bki %
tfs1 = 0:fs1:T; 
tfs2 = 0:fs2:T;
tfs3 = 0:fs3:T;

% Amplitudy sinusoid %
A1 = 50;
A2 = 100;
A3 = 150;

% Tworzenie sygna?u z sumy sinus?w %
s1 = A1 * sin(2*pi*f1*tfs1);
s2 = A2 * sin(2*pi*f2*tfs2);
s3 = A3 * sin(2*pi*f3*tfs3);

x = s1 + s2 + s3; % sygna? x z sumy sinus?w %

figure(1);
subplot(2,1,1);
hold all;
plot(s1, 'r-o');
plot(s2, 'g-o');
plot(s3, 'b-o');
title('Trzy sinusy do sumowania');
legend('50Hz','100Hz','150Hz');

figure(1);
subplot(2,1,2);
plot(x, 'r-x')
title('Zsumowane sinusy');
legend('50Hz + 100Hz + 150Hz');

% Budujemy macierze DCT i IDCT %
N = 101;
s = sqrt(1/N);
for k = 1:N
    for n = 1:N
        A(k,n) = s * cos(pi*(k-1)/N *((n-1)+0.5));
    end
    s = sqrt(2/N);
end
%IDCT%
S = transpose(A);

% Wy?wietlanie wierszy A i kolumn S %

%{
figure(2);
for i=1:100
    plot(A(i,:),'r-x');
    hold on;
    plot(S(:,i),'b-o');
    title('Wiersze i kolumny DCT/IDCT');
    legend('wiersze DCT','kolumny IDCT');
    fprintf('WIERSZ/KOLUMNA NR: %u\n', i);
    pause
    hold off;
end
%}
% Analiza y=Ax %
y = A*x'; 
rekonstr = S*y;

display('Warto?ciy(1:N)');
display(y(1:N));

% usuwanie wspolczynnikow zerowych %
pom = 1;
for p = 1:101
    if abs(y(p))>(10^(-11));
        ynz(pom)=y(p);
        pom = pom+1;
    end
end




% ANALIZA CZ?STOTLIWO?CI WYKRESY %
% por?wnaj warto?ci wsp?czynnik?w niezerowych z warto?ciami amplitud %
% sk?adowych sygna?u %
f =(0:N-52)*fs/N;
figure(3);
subplot(3,2,1);
stem(f,ynz, 'b');
title('Obserwacja wsp?czynnik?w');
xlabel('Cz?stotliwo?? [Hz]');
subplot(3,2,2);
plot(x, 'r-o');
title('Sygna? sumy sinus?w');

% ZMIANA f2 na 105Hz%
f2 = 105;

s2 = A2 * sin(2*pi*f2*tfs2);
x2 = s1 + s2 + s3; % sygna? x z sumy sinus?w %
y = A*x2'; 
rekonstr2 = S*y;



f =(0:N-1)*500/N;
subplot(3,2,3);
stem(f,y, 'b');
title('Obserwacja wsp?czynnik?w f2 = 105Hz');
xlabel('Cz?stotliwo?? [Hz]');
subplot(3,2,4);
plot(x2, 'r-o');
title('Sygna? sumy sinus?w f2 = 105Hz');

% ZMIANA wszystkich f o 2.5Hz %
f1 = 52.5; 
f2 = 102.5;
f3 = 152.5;

s1 = A1 * sin(2*pi*f1*tfs1);
s2 = A2 * sin(2*pi*f2*tfs2);
s3 = A3 * sin(2*pi*f3*tfs3);
x3 = s1 + s2 + s3; % sygna? x z sumy sinus?w %
y = A*x3'; 
rekonstr3 = S*y;

f =(0:N-1)*500/N;
subplot(3,2,5);
stem(f,y, 'b');
title('Obserwacja wsp?czynnik?w f + 2.5Hz');
xlabel('Cz?stotliwo?? [Hz]');
subplot(3,2,6);
plot(x3, 'r-o');
title('Sygna? sumy sinus?w f + 2.5Hz');

% POR?WNANIE ODTWORZE? SYGNA??W %
figure(4);
subplot(3,2,1);
plot(x, 'b');
title('Sygna? suma sinus?w');
subplot(3,2,2);
plot(rekonstr, 'b');
title('Odtworzony sygna? sumy sinus?w');

subplot(3,2,3);
plot(x2, 'b');
title('Sygna? suma sinus?w f2=105Hz');
subplot(3,2,4);
plot(rekonstr2, 'b');
title('Odtworzony sygna? sumy sinus?w f2=105Hz');

subplot(3,2,5);
plot(x3, 'b');
title('Sygna? suma sinus?w f + 2.5Hz');
subplot(3,2,6);
plot(rekonstr3, 'b');
title('Odtworzony sygna? sumy sinus?w f + 2.5Hz');