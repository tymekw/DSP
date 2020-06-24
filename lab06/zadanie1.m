clc; clear all; close all;
load('butter.mat');
fs=16000;
f = 0:1:fs/2-1;
w =2*pi*f;

[zz,pp,kk] = bilinear(z,p,k,fs);
f1 = 1189;
f2 = 1229;
bw = f2-f1; %szerokosc pasma
%bA = k* poly(z);
%aA = poly(p);
[bA,aA]=zp2tf(z,p,k);%Convert zero-pole-gain filter parameters to transfer function form
HA = freqs(bA,aA,w);

%bC = kk* poly(zz);
%aC = poly(pp);
[bC,aC] = zp2tf(zz,pp,kk);
HC = freqz(bC,aC,f,fs);   
save('f_cyfrowy.mat','HC');
figure;
plot(f,20*log10(abs(HC)),'r');
hold on;
plot(f,20*log10(abs(HA)),'b');
xline(f1, 'g');
xline(f2, 'g');
yline(-3, 'g');
xlim([1160, 1250]);
ylim([-6,0]);
title("|H(w)| [dB]")
legend("cyfrowy", "analogowy");
grid;


% druga czêœæ -> filtracja WYK£AD


fs=16000;
T = 1; %czas trwania sygna³u
Nx = T*fs; %ilosc probek
dt =1/fs; t = dt*(0:Nx-1); %odstêp miedzy i do sinusa
x = sin(2*pi*1209*t)+sin(2*pi*1272*t);
b = bC;
a = aC;

M = length(b);
N = length(a); a=a(2:N); N=N-1;

bx = zeros(1,M); %bufor x
by = zeros(1,N); %bufor y

for n =1 :Nx
    bx = [x(n) bx(1:M-1)]; %aktualizacja bufora; obecna probka na poczatek
    y(n) = sum(bx .* b) - sum(by .*a);
    by = [y(n) by(1:N-1)]; %y do bufora 
end


yy = filter(bC,aC,x);

figure
subplot(311); plot(t,x); title('x(n)'); grid;
subplot(312); plot(t,y); title('y(n)'); grid;
subplot(313); plot(t,yy); title('y(n) filter()'); grid;
%Wykres sygnalow w czasie
if(0)
figure(3)
plot(t,x,'b-x');
hold on;
plot(t,real(y),'r-o');
hold on;
plot(t,yy,'m-*');
title('Sygnal:');
xlabel('Czas [s]');
ylabel('Amplituda');
legend('orginalny sygnal','filtrowany cyfrowo','filter()');
xlim([0.5 0.515])
end

figure
%pozbywam siê stanów przejœciowych 
k = Nx/2+1:Nx; f0 = fs/(Nx/2); f = f0*(0:Nx/2-1); 
subplot(311); plot(f,20*log10(abs(2*fft(x(k)))/(Nx/2))); title('x(n)'); grid;
subplot(312); plot(f,20*log10(abs(2*fft(y(k)))/(Nx/2))); title('y(n)'); grid;
subplot(313); plot(f,20*log10(abs(2*fft(yy(k)))/(Nx/2))); title('y(n) filter()'); grid;

%Wykres widma sygnalow
figure
plot(f,20*log10(abs(2*fft(x(k)))/(Nx/2)),'k-x');
hold on;
plot(f,20*log10(abs(2*fft(y(k)))/(Nx/2)),'r-*');
hold on;
plot(f,20*log10(abs(2*fft(yy(k)))/(Nx/2)), 'b-o')
title('Widma sygnalu:');
xlabel('Czestotliwosc [Hz]');
ylabel('Tlumienie [dB]');
legend('orginalny sygnal','filtrowany cyfrowo','filter()');
xlim([1150 1300]);
