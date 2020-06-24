clear all;
close all;

%% za�adowanie sygna��w
load('lab08_am.mat');
x=s0;

%% wygenerowanie teoretycznej odpowiedzi impulsowej
fs=1000;        %czestotliwosc probkowania
fc=200;         %czestotliwosc nosna
M=64;           %polowa dlugosci filtra
N=2*M+1;
n=1:M;
h=(2/pi)*sin(pi*n/2).^2 ./n;    %po�owa odpowiedzi impulsowej (TZ str. 352)
h=[-h(M:-1:1) 0 h(1:M)];        %ca�a odpowied� dla n = ?M,...,0,...,M

%% wymna�amy przez okno Blackmana
w=blackman(N); w=w';            
hw=h.*w; % wymno�enie odpowiedzi impulsowej z oknem

%% widmo Fouriera oraz wykresy
%{
m = -M : 1 : M; % dla filtra nieprzyczynowego (bez przesuni�cia o M pr�bek w prawo)
% m = 0 : N-1; % dla filtra przyczynowego (z przesuni�ciem o M pr�bek w prawo)
NF=500; fn=0.5*(1:NF-1)/NF;
for k=1:NF-1
H(k)=sum( h .* exp(-j*2*pi*fn(k)*m) );
HW(k)=sum( hw .* exp(-j*2*pi*fn(k)*m) );
end

figure(1);
subplot(2,2,1);
stem(m,h); grid; title('h(n)'); xlabel('n');
subplot(2,2,2);
stem(m,hw); grid; title('hw(n)'); xlabel('n'); 
subplot(2,2,3);
plot(fn,abs(H)); grid; title('|H(fn)|'); xlabel('f norm]'); 
subplot(2,2,4);
plot(fn,abs(HW)); grid; title('|HW(fn)|'); xlabel('f norm]');
%}

%% por�wnanie z funkcja hilbert() matlaba
HIL = hilbert(x);

%% filtracja odpowiedzia imp.
y=conv(x,hw); % filtracja sygna�u x(n) za pomoc� odp. impulsowej hw(n); otrzymujemy Nx+N?1 pr�bek
y=y(N:1000);  % odci�cie stan�w przej�ciowych (po N?1 pr�bek) z przodu sygna�u y(n)
c=x(M+1:1000-M);   % odci�cie tych pr�bek z x(n), dla kt�rych nie ma poprawnych odpowiednik�w w y(n)
m=sqrt(c.^2+y.^2);  %Obwiednia to pierwiastek z sumy kwadrat�w sygna��w x i jego transformacji Hilberta HT(x).
NFFT=2^nextpow2(fs);
Y=fft(m,NFFT)/fs; %transformata fouriera syg m
f=fs/2*linspace(0,1,NFFT/2+1);
figure(2);
plot(f,2*abs(Y(1:NFFT/2+1)));
title('FFT obwiedni');

%% odczytane parametry sygna�u moduluj�cego
f1=2;
f2=40;
f3=60;
A1=0.3355;
A2=0.1738;
A3=0.2032;

%% sygna� moduluj�cy - suma 3 cosinusow
t=0:1/fs:1-1/fs;
xr=1+A1*cos(2*pi*f1*t)+A2*cos(2*pi*f2*t)+A3*cos(2*pi*f3*t); %sygna� moduluj�cy
xr=xr(M+1:1000-M); %odci�cie pr�bek

%% sygna� z pliku po transformacie hilberta vs sygna� skonstruowany
figure(3);
hold all;
plot(m, 'b');
plot(xr,'k');
plot(abs(HIL),'y');
title('Por�wnanie sygna��w');
legend('FIR - nasz hilbert','odtworzone suma cos','matlab hilbert()'); %dla przejrzystosci mozna wywalic hilberta matlaba

%% sygna� z pliku zmodulowany vs sygna� skonstruowany zmodulowany
figure(4)
xp=sin(2*pi*fc*t);  %no�na z tre�ci
xp=xp(M+1:1000-M);  %odci�cie pr�bek
xz=xp.*xr;          %sygna� zmodulowany po odtworzeniu
hold on;
plot(xz,'b');
plot(y,'r');
legend('Sygna� skonstruowany i zmod.','Sygna� zmod. z pliku');
