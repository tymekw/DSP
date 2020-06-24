clc;
clear all;
close all;
[s, fs] = audioread('mowa_3.wav'); s=s';
[sA,fs] = audioread('mowa_1.wav');   sA=sA';
[sB,fs] = audioread('mowa_2.wav');   sB=sB';

x = sB; %signal echo
d = s; %disturbed signal
s = sA; %desired signal
Nx = length(x);
f=0:fs/500:fs/2;
dt = 1/fs;
t = dt*(0:Nx-1);


M=10;
mi=0.83;

% Calculation of optimal Wiener filter and limits for mi
for k = 0 : M
  rxx(k+1) = sum( x(1+k:Nx) .* x(1:Nx-k) )/(Nx-M);  % auto-correlation of x(n)
  rdx(k+1) = sum( d(1+k:Nx) .* x(1:Nx-k) )/(Nx-M);  % cross-correlation of d(n) and x(n)
end
Rxx = toeplitz(rxx,rxx);         % symmetrical autocorrelation matrix of x(n)
h_opt = Rxx\rdx';                % weights of the optimal Wiener filter
lambda = eig( Rxx );             % eigenvalue decomposition of Rxx
lambda = sort(lambda,'descend'); % sorting eigenvalues
%disp('Suggested values of mi')
%mi1_risc = 1/lambda(1),
%figure;
%subplot(211); stem( h_opt ); grid; title('Optimal Wiener filter h(n)');
%subplot(212); stem( lambda ); grid; title('Eigenvalues of matrix Rxx');



%adaptive filtration
bx = zeros(M,1); %buffor
h = zeros(M,1); %filter weights
%h = h_opt(1:15);

y=[];
e=[];
Rinv = 0.5*eye(M,M); 
for n=1:length(x)
    bx = [x(n);bx(1:M-1)];
    y(n) = h'*bx; %FIR
    e(n) = d(n)-y(n); %liczymy e
    h = h+mi*e(n)*bx; %LMS
    if(0)
        lambd = 0.999; 
        Rinv = (Rinv - Rinv*bx*bx'*Rinv/(lambd+bx'*Rinv*bx))/lambd;
        h = h + (e(n) * Rinv * bx );
    end
end


soundsc(e,fs)
figure; 
plot(t,s, 'r*-',t,e,'b'); grid; 
xlabel('czas [s]');
title('orginalny, po filtracji');legend("orginalny","mój")

figure; 
stem(1:M+1,h_opt); hold on;
stem(1:M,h); grid;
xlabel('n');
title('wagi filtrów'); legend("optymalny","mój");

figure; 
plot(f,20*log10(abs(freqz(h_opt(1:M),1,f,fs)))); hold on;
plot(f,20*log10(abs(freqz(h,1,f,fs)))); grid;
title('odpowiedz czestotliwosciowa'); legend("optymalny","mój");
xlabel('Hz');


f = (0:1:Nx-1)*fs/Nx
figure;
plot(f,20*log10(abs(fft(d))));
hold all;
plot(f,20*log10(abs(fft(s))));
plot(f,20*log10(abs(fft(e))));
legend("zaszumiony","oryginalny","odszumiony");
title("fft(AWGN 20dB)");


SNR = 10*log10((1/Nx * sum(s.^2))/(1/Nx*sum((s-e).^2)));









