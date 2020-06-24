clc;
clear all; close all;

fs = 8000;
T=1;
dt = 1/fs;
t = (0:dt:T-dt);
A1=-0.5; A2=1;
f1=34.2; f2=115.5;
f = 0:fs/500:fs/2;



dref = A1*cos(2*pi*t*f1)+A2*cos(2*pi*t*f2); %sygna³ czysty

%moc drefu-> suma kwadratow/ fs
drefpwr = 10*log10(sum(dref.^2)/fs);

d10 = awgn(dref,10,'measured'); %WE: sygna³ odniesienie dla x
d20 = awgn(dref,20,'measured'); %WE: sygna³ odniesienie dla x
d40 = awgn(dref,40,'measured'); %WE: sygna³ odniesienie dla x

x10 = [d10(1) d10(1:end-1)];   %WE sygna³ filtrowany, opozniony d
x20 = [d20(1) d20(1:end-1)]; 
x40 = [d40(1) d40(1:end-1)]; 

M = 89; %dlugosc filtra
mi =0.009; %0.047; % wspolczynnik szybkosci adaptacji


y10=[]; e10=[]; %sygnaly na wyjsciu
bx10 = zeros(M,1); %bufor na probki wejsciowe x
h10 = zeros(M,1); %poczatkowe puste wagi
figure;
Rinv = 0.5*eye(M,M);
for n=1:length(x10)
    bx10 = [x10(n);bx10(1:M-1)]; %nowa probka do bufora
    y10(n) = h10'*bx10; %FIR
    e10(n) = d10(n)-y10(n); %liczymy e
    h10 = h10+mi*e10(n)*bx10; %LMS
    %h10 = h10+mi*e10(n) *bx10/(bx10'*bx10); %NLMS
    if(0)
        lambd = 0.999; 
        Rinv = (Rinv - Rinv*bx10*bx10'*Rinv/(lambd+bx10'*Rinv*bx10))/lambd;
        h10 = h10 + (e10(n) * Rinv * bx10 );
    end
    if(0) 
        subplot(211); stem(h10); xlabel('n'); title('h10(n)'); grid;
        subplot(212); plot(f,abs(freqz(h10,1,f,fs))); xlabel('f (Hz)');
        title('|H(f)|'); grid; % pause
        pause;
    end
      
end

Rinv = 0.5*eye(M,M);
y20=[]; e20=[]; %sygnaly na wyjsciu
bx20 = zeros(M,1); %bufor na probki wejsciowe x
h20 = zeros(M,1); %poczatkowe puste wagi
for n=1:length(x20)
    bx20 = [x20(n);bx20(1:M-1)]; %nowa probka do bufora
    y20(n) = h20'*bx20; %FIR
    e20(n) = d20(n)-y20(n); %liczymy e
    h20 = h20+mi*e20(n)*bx20; %LMS
    %h20 = h20+mi*e20(n) *bx20/(bx20'*bx20); %NLMS
    if(0)
        lambd = 0.999; 
        Rinv = (Rinv - Rinv*bx20*bx20'*Rinv/(lambd+bx20'*Rinv*bx20))/lambd;
        h20 = h20 + (e20(n) * Rinv * bx20 );
    end
end


Rinv = 0.5*eye(M,M);
y40=[]; e40=[]; %sygnaly na wyjsciu
bx40 = zeros(M,1); %bufor na probki wejsciowe x
h40 = zeros(M,1); %poczatkowe puste wagi
for n=1:length(x40)
    bx40 = [x40(n);bx40(1:M-1)]; %nowa probka do bufora
    y40(n) = h40'*bx40; %FIR
    e40(n) = d40(n)-y40(n); %liczymy e
    h40 = h40+mi*e40(n)*bx40; %LMS
    %h40 = h40+mi*e40(n) *bx40/(bx40'*bx40); %NLMS
    if(0) %RLS
        lambd = 0.999; 
        Rinv = (Rinv - Rinv*bx40*bx40'*Rinv/(lambd+bx40'*Rinv*bx40))/lambd;
        h40 = h40 + (e40(n) * Rinv * bx40 );
    end
end

figure;
plot(t,dref);
hold all;
plot(t,d10);
plot(t,y10);
legend("oryginalny","zaszumiony","odszumiony");
title("AWGN 10dB");

figure;
plot(t,dref);
hold all;
plot(t,d20);
plot(t,y20);
legend("oryginalny","zaszumiony","odszumiony");
title("AWGN 20dB");

figure;
plot(t,dref);
hold all;
plot(t,d40);
plot(t,y40);
legend("oryginalny","zaszumiony","odszumiony");
title("AWGN 40dB");

licznik=0;
mianownik10 =0;
mianownik20=0;
mianownik40=0;
for n = 200:length(dref) %%obcinam dla SNR !!!!!!!!!!!!!
    licznik = licznik + dref(n)^2;
    mianownik10 = mianownik10 + (dref(n) - y10(n))^2;
    mianownik20 = mianownik20 + (dref(n) - y20(n))^2;
    mianownik40 = mianownik40 + (dref(n) - y40(n))^2;
end

SNR10 = 10*log10(((1/length(d10))*licznik)/((1/length(d10))*mianownik10));
SNR20 = 10*log10(((1/length(d20))*licznik)/((1/length(d20))*mianownik20));
SNR40 = 10*log10(((1/(length(d40)-200))*licznik)/((1/(length(d40)-200))*mianownik40));

% Calculation of optimal Wiener filter and limits for mi for d10
Nx = 1*fs;
for k = 0 : M
  rxx(k+1) = sum( x10(1+k:Nx) .* x10(1:Nx-k) )/(Nx-M);  % auto-correlation of x(n)
  rdx(k+1) = sum( d10(1+k:Nx) .* x10(1:Nx-k) )/(Nx-M);  % cross-correlation of d(n) and x(n)
end
Rxx = toeplitz(rxx,rxx);         % symmetrical autocorrelation matrix of x(n)
h_opt = Rxx\rdx';                % weights of the optimal Wiener filter
lambda = eig( Rxx );             % eigenvalue decomposition of Rxx
lambda = sort(lambda,'descend'); % sorting eigenvalues
disp('Suggested values of mi')
mi_risc = 1/sum(lambda),

figure;
stem(1:M+1, h_opt); hold on;
stem(1:M,h10); grid;
title("wagi filtrów");legend("optymalny","mój");xlabel("n");
%fvtool(h10)
figure; subplot(111); plot(f,20*log10(abs(freqz(h_opt(1:M),1,f,fs))),f,20*log10(abs(freqz(h10,1,f,fs)))); grid;
title('odpowiedz czestotliwosciowa');legend("optymalny","mój");   xlabel('Hz');

f = 0:1:fs-1
figure;
plot(f,20*log10(abs(fft(dref))));
hold all;
plot(f,20*log10(abs(fft(d20))));
plot(f,20*log10(abs(fft(y20))));
legend("oryginalny","zaszumiony","odszumiony");
title("fft(AWGN 20dB)");

