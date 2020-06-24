close all;
clear all;
load('butter.mat'); % z, p, k
 
%Dane
% fp1 = 1189; % czestotliwosc graniczna dolna
% fp2 = 1229; % czestotliwosc graniczna gorna
N = 16000;  % ilosc probek gdy 1s i 16khz
fs = 16000; % czestotliwosc probkowania
 
A1 = 1; % amplitudy do generowanego sygnalu
A2 = 1;
f1 = 1209; % czestotliwosci do generowanego sygnalu
f2 = 1272;
 
f = linspace(0,N,N); %0:1:16000
w = f*2*pi; %omega
 
 
%Charakterystyka filtru w zadaniu
%Convert zero-pole-gain filter parameters to transfer function form
[b_butter_zadany, a_butter_zadany] = zp2tf(z, p, k);
%odpowiedz impulsowa
[H_butter, W_butter] = freqs(b_butter_zadany, a_butter_zadany, w);
 
%Transformata biliniowa (konwersja H(s)->H(z))
[b_trans_bilin,a_trans_bilin] = bilinear(b_butter_zadany,a_butter_zadany,fs);
[H_po_bilin,W_po_bilin] = freqz(b_trans_bilin,a_trans_bilin,f,fs);
save('f_cyfrowy.mat','b_trans_bilin', "a_trans_bilin");
 
%Zamiana na wartosci logarytmiczne
H_butter_abs = abs(H_butter);
H_bilin_abs = abs(H_po_bilin);
H_butter_log = 20*log10(H_butter_abs);
H_bilin_log = 20*log10(H_bilin_abs);
 
%Wykres transmitancji filtrow w f
x3dB = [1000 1500];
y3dB = [-3 -3];
 
figure(1);
hold all;
plot(f,H_butter_log,'ko-'); %1168-1205.8
plot(f,H_bilin_log,'b*-'); %1188.9-1228.9 %przesuniecie przez nielinilwosc
plot(x3dB,y3dB,'r');
title('Transmitancja filtrów:');
xlabel('Czestotliwosc [Hz]');
ylabel('Tlumienie [dB]');
legend('analogowy','cyfrowy');
xlim([1160 1270]);
ylim([-6 0]);
% WYKRES PE£NEJ TRANSMITANCJI
% hold all;
% plot(f,H_butter_log,'ko-'); %1168-1205.8
% plot(f,H_bilin_log,'b*-'); %1188.9-1228.9
% plot(x3dB,y3dB,'r');
% title('Transmitancja filtrów:');
% xlabel('Czestotliwosc [Hz]');
% ylabel('Tlumienie [dB]');
% legend('analogowy','cyfrowy');
 
%Sygnal i jego filtracja
t = 1/fs:1/fs:1;
x_sygnal = A1*sin(2*pi*f1*t) + A2*sin(2*pi*f2*t);
 
X_sygnal = fft(x_sygnal)/(N/2); % przejscie na dziedzine czestotliwosc
X_sygnal_abs = abs(X_sygnal);
X_sygnal_log = 20*log10(X_sygnal_abs);
 
Y_sygnal = H_po_bilin .* X_sygnal;
Y_sygnal_abs = abs(Y_sygnal);
Y_sygnal_log = 20*log10(Y_sygnal_abs);
y_sygnal = ifft(Y_sygnal); % powrot do czasu
 
% porownanie
y_filter = filter(b_trans_bilin,a_trans_bilin,x_sygnal);
Y_filter = 20*log10(abs(fft(y_filter)/(N/2)));
 
%Wykres widma sygnalow
figure(2)
plot(f,X_sygnal_log,'k-x');
hold on;
plot(f,Y_sygnal_log,'r-*');
hold on;
plot(f, Y_filter, 'b-o')
title('Widma sygnalu:');
xlabel('Czestotliwosc [Hz]');
ylabel('Tlumienie [dB]');
legend('orginalny sygnal','filtrowany cyfrowo','filter()');
xlim([1160 1270]);
 
%Wykres sygnalow w czasie
figure(3)
plot(t,x_sygnal,'b-x');
hold on;
plot(t,N*real(y_sygnal),'r-*');
hold on;
plot(t,y_filter,'m-*');
title('Sygnal:');
xlabel('Czas [s]');
ylabel('Amplituda');
legend('orginalny sygnal','filtrowany cyfrowo','filter()');
xlim([0.5 0.515])