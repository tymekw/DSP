clc;
clear all; close all;
%BP = LP2 - LP1

fpr = 1200; %Hz
df = 200; %Hz pasmo przepustowe
fc = 300; %Hz czestotliwosc srodkowa pasma przepustowego
f1 = 300-(200/2); %Hz lewa strona
f2 = 300+(200/2); %Hz prawa strona
%MM = 128; %d³ugoœæ flitru z³a bo parzysta
MM = 129; %d³ugoœæ filtru poprawna
M = (MM-1)/2;
n = -M:1:M; % filt przechodz¹cy przez zero

%LP filter do 400Hz
hLP2 = 2*f2/fpr*sin(2*pi*f2/fpr*n)./(2*pi*f2/fpr*n);
hLP2(M+1) = 2*f2/fpr;

%LP filter do 200Hz
hLP1 = 2*f1/fpr*sin(2*pi*f1/fpr*n)./(2*pi*f1/fpr*n);
hLP1(M+1) = 2*f1/fpr;

%BP filter
h = hLP2-hLP1;

%BP + okna
h1 = h .* boxcar(MM)';
h2 = h .* hanning(MM)';
h3 = h .* hamming(MM)';
h4 = h .* blackman(MM)';
h5 = h .* blackmanharris(MM)';


%odpowiedŸ impulsowa filtrów
f=0:1:fpr/2;
H1 = freqz(h1,1,f,fpr);   %H = polywal(h,z) z = roots(h)
H2 = freqz(h2,1,f,fpr);
H3 = freqz(h3,1,f,fpr);
H4 = freqz(h4,1,f,fpr);
H5 = freqz(h5,1,f,fpr);

%char a-cz filtru
figure;
plot(f,20*log10(abs(H1)));
hold on;
plot(f,20*log10(abs(H2)));
plot(f,20*log10(abs(H3)));
plot(f,20*log10(abs(H4)));
plot(f,20*log10(abs(H5)));
xlabel('f [Hz]'); title("|H(f)|"); grid;
legend("Prostok¹tny", "Hanning","Hamming","Blackman","Blackman-Harris");

%char fa-cz filtru
figure;
hold on;
plot(f,unwrap(angle(H1)));
plot(f,unwrap(angle(H2)));
plot(f,unwrap(angle(H3)));
plot(f,unwrap(angle(H4)));
plot(f,unwrap(angle(H5)));
xlabel('f [Hz]'); title("Phase H(f)"); grid;
legend("Prostok¹tny", "Hanning","Hamming","Blackman","Blackman-Harris");


%budowa sygna³u
T = 1;
Nx = T*fpr;
dt =1/fpr; t = dt*(0:Nx-1); %odstêp miedzy i do sinusa
x = sin(2*pi*300*t)+sin(2*pi*500*t)+sin(2*pi*100*t);

%filtracja y = conv(x,h); y = filter(h,1,x);
Nh = length(h);
%Prostok¹tny
hx1 = zeros(1,Nh);
for k = 1:Nx
    hx1 = [x(k) hx1(1:Nh-1)];
    y1(k) = sum(hx1 .* h1);
end

%Hanning
hx2 = zeros(1,Nh);
for k = 1:Nx
    hx2 = [x(k) hx2(1:Nh-1)];
    y2(k) = sum(hx2 .* h2);
end

%Hamming
hx3 = zeros(1,Nh);
for k = 1:Nx
    hx3 = [x(k) hx3(1:Nh-1)];
    y3(k) = sum(hx3 .* h3);
end

%Blackman
hx4 = zeros(1,Nh);
for k = 1:Nx
    hx4 = [x(k) hx4(1:Nh-1)];
    y4(k) = sum(hx4 .* h4);
end

%Blackman-Harris
hx5 = zeros(1,Nh);
for k = 1:Nx
    hx5 = [x(k) hx5(1:Nh-1)];
    y5(k) = sum(hx5 .* h5);
end

%figure;
%subplot(211); plot(t,x); title('x(n)'); grid;
%subplot(212); plot(t,y1); title('y1(n)'); grid;


%tylko prawa strona bo to kopia lewej
k=Nx/2+1:Nx;
f=(fpr/(Nx/2))*(0:Nx/2-1);

figure;
plot(f,20*log10(2*abs(fft(x(k)))/(Nx/2)));
title('|X(n)|'); grid;

figure;
hold on;
plot(f,20*log10(2*abs(fft(y1(k)))/(Nx/2)));
plot(f,20*log10(2*abs(fft(y2(k)))/(Nx/2)));
plot(f,20*log10(2*abs(fft(y3(k)))/(Nx/2)));
plot(f,20*log10(2*abs(fft(y4(k)))/(Nx/2)));
plot(f,20*log10(2*abs(fft(y5(k)))/(Nx/2)));
legend("Prostok¹tny", "Hanning","Hamming","Blackman","Blackman-Harris");
title('|Y(n)|'); grid;


figure;
psd(spectrum.welch("Hamming",1024),y1,'Fs',fpr);
title('prostokatny'); grid;

figure;
psd(spectrum.welch("Hamming",1024),y2,'Fs',fpr);
title('Hanning'); grid;

figure;
psd(spectrum.welch("Hamming",1024),y3,'Fs',fpr);
title('hamming'); grid;

figure;
psd(spectrum.welch("Hamming",1024),y4,'Fs',fpr);
title('blackman'); grid;

figure;
psd(spectrum.welch("Hamming",1024),y5,'Fs',fpr);
title('blackman-harris'); grid;
%legend("Prostok¹tny", "Hanning","Hamming","Blackman","Blackman-Harris");

figure;
psd(spectrum.welch("Hamming",1024),x,'Fs',fpr);
title('wejsciowy'); grid;

