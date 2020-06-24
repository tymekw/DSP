clc;
clear all;
N = 100;
for k = 0:N-1
    for n = 0:N-1
        A(k+1,n+1) = (1/N)*((exp(1i*2*pi/N))^(-k*n));
    end
end

fs = 1000;   %Hz
f1 = 100;    %Hz
f2 = 125;    %Hz
A1 = 1;    
A2 = 0.0001;
fi1 = pi/7;
fi2 = pi/11;

dt = 1/fs;  %okresy próbkowañ
t = (0:N-1)*dt; %momenty próbek

%sygna³
x1 = A1*cos(2*pi*f1*t + fi1);
x2 = A2*cos(2*pi*f2*t + fi2);
x = x1+x2;

%DtFT 
f = 0:0.1:500;
for k = 1:length(f)
    X(k) = 1/N*(sum(x.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
end

%fx = fs*(0:50*N)/100*N;

figure;
plot(f,abs(X));
title("abs(X)");


% prostok¹tnym, Hamminga, Blackmana, Czebyszewa (t³umienie 100
%dB) i Czebyszewa (t³umienie 120 dB)

%prostok¹tny
x1 = x.* boxcar(N)';
x2 = x.* hamming(N)';
x3 = x.* blackman(N)';
x4 = x.* chebwin(N,100)';
x5 = x.* chebwin(N,120)';

%DtFT dla okien
for k = 1:length(f)
    X1(k) = 1/N*(sum(x1.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
    X2(k) = 1/N*(sum(x2.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
    X3(k) = 1/N*(sum(x3.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
    X4(k) = 1/N*(sum(x4.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
    X5(k) = 1/N*(sum(x5.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
end

%wykresy widm
figure;
hold on;

plot(f,20*log10(abs(X1)));
plot(f,20*log10(abs(X2)));
plot(f,20*log10(abs(X3)));
plot(f,20*log10(abs(X4)));
plot(f,20*log10(abs(X5)));
legend("prostok¹t","Hamming","Blackman","Czebyszew 100", "Czebyszew 120");

%okna:
x10 = boxcar(N)';
x20 = hamming(N)';
x30 = blackman(N)';
x40 = chebwin(N,100)';
x50 = chebwin(N,120)';

figure;
hold on;
plot(x10);
plot(x20);
plot(x30);
plot(x40);
plot(x50);
legend("prostok¹t","Hamming","Blackman","Czebyszew 100", "Czebyszew 120");


%=========================================================================%
N=1000;
for k = 0:N-1
    for n = 0:N-1
        A(k+1,n+1) = (1/N)*((exp(1i*2*pi/N))^(-k*n));
    end
end

fs = 1000;   %Hz
f1 = 100;    %Hz
f2 = 125;    %Hz
A1 = 1;    
A2 = 0.0001;
fi1 = pi/7;
fi2 = pi/11;

dt = 1/fs;  %okresy próbkowañ
t = (0:N-1)*dt; %momenty próbek

%sygna³
x1 = A1*cos(2*pi*f1*t + fi1);
x2 = A2*cos(2*pi*f2*t + fi2);
x = x1+x2;

%DtFT 
f = 0:0.1:500;
for k = 1:length(f)
    X(k) = 1/N*(sum(x.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
end

%fx = fs*(0:5*N)/10*N;

figure;
plot(f,abs(X));



% prostok¹tnym, Hamminga, Blackmana, Czebyszewa (t³umienie 100
%dB) i Czebyszewa (t³umienie 120 dB)

%prostok¹tny
x1 = x.* boxcar(N)';
x2 = x.* hamming(N)';
x3 = x.* blackman(N)';
x4 = x.* chebwin(N,100)';
x5 = x.* chebwin(N,120)';

%DtFT dla okien
for k = 1:length(f)
    %X1(k) = 1/N*(sum(x1.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
    %X2(k) = 1/N*(sum(x2.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
    %X3(k) = 1/N*(sum(x3.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
    X4(k) = 1/N*(sum(x4.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
    X5(k) = 1/N*(sum(x5.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
end

%wykresy widm
figure;
hold on;

%plot(fx,abs(X1));
%plot(fx,abs(X2));
%plot(fx,abs(X3));
plot(f,20*log10(abs(X4)));
plot(f,20*log10(abs(X5)));
legend("Czebyszew 100", "Czebyszew 120");



