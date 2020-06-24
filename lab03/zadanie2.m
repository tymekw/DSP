clc;
clear all;

N=100;
for k = 0:N-1
    for n = 0:N-1
        A(k+1,n+1) = (1/N)*((exp((1i*2*pi)/N))^(-k*n)); %dziele przez N tu
    end
end


fs = 1000;   %Hz
f1 = 125;    %Hz
f2 = 200;    %Hz
A1 = 100;    
A2 = 200;
fi1 = pi/7;
fi2 = pi/11;

dt = 1/fs;  %okresy próbkowañ
t = (0:N-1)*dt; %momenty próbek

%sygna³
x1 = A1*cos(2*pi*f1*t + fi1);
x2 = A2*cos(2*pi*f2*t + fi2);
x = x1+x2;

%DFT (X = Ax)
X = A*x';


%=========================================================================%

%DFT
X1=X;

%dodanie M zer
M = 100;
xz = [x,zeros(1,M)];

% skalowanie
X2 = fft(xz)./(N+M);

%oblicznie X3 ze wzoru na DtFT
f = 0:0.25:1000;
for k = 1:length(f)
    X3(k) = 1/N*(sum(x.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
end


%DFT
fx1 = fs*(0:N-1)/N;
fx2 = fs*(0:M+N-1)/(N+M);
fx3 = fs*(0:40*N)/(40*N);

figure;
plot(fx1,X,'o',fx2,X2,'bx',f,X3,'k-');
legend("DFT","DFT z zerami","DtFT")



%=========================================================================%

%oblicznie X3 ze wzoru na DtFT dla innego f
f=-2000:0.25:2000;
for k = 1:length(f)
    X3(k) = 1/N*(sum(x.*exp(-1i*2*pi*(f(k)/fs) * (0:N-1))));
end

%DFT
fx1 = fs*(0:N-1)/N;
fx2 = fs*(0:M + N-1)/(N+M);
fx3 = fs*(0:160*N)/(160*N);

figure;
plot(fx1,X,'o',fx2,X2,'bx',f,X3,'k-');
legend("DFT","DFT z zerami","DtFT")







