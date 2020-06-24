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

%rysowanie widma
rex = real(X);
imx = imag(X);
modx = abs(X);
phx = atan(imx./rex);


f = fs*(0:N-1)/N;

figure;
subplot(411);
stem(f,rex);
title("rzeczywista");
subplot(412);
stem(f,imx);
title("urojona");
subplot(413);
stem(f,modx);
title("modu³");
subplot(414);
stem(f,phx);
title("faza");

%=========================================================================%

%macierz rekonstrukcji -> zamienione rzêdy na kolumny
B = A';
xr = B*X;
xr = xr';

%czy x == xr?
maxBlad = max(abs(x-xr));
fprintf("maksymalny blad rekonstrukcji sygna³u DFT:  %d \n" ,maxBlad);

%rekonstrukcja fft
X1 = fft(x);
xr1 = ifft(X1);
maxBlad1 = max(abs(x-xr1));
fprintf("maksymalny blad rekonstrukcji sygna³u FFT:  %d \n" ,maxBlad1);

%roznica
roznica = abs(x-xr1);
fprintf("max roznica miedzy metodami: %d \n" ,max(roznica));

%figure;
%plot(roznica);
%=========================================================================%
%in DFT f do porownania te w macierzy A
%musz¹ byæ k-ta wielokrotnoœci¹ fs/N k=0:1:N-1
% fs/N = 10 f1 nie wystêpuje najblizej jest 120 i 130 


f1 = 125;    %Hz

%sygna³
x1 = A1*cos(2*pi*f1*t + fi1);
x2 = A2*cos(2*pi*f2*t + fi2);
x = x1+x2;

%DFT (X = Ax)
X = A*x';

%rysowanie widma
rex = real(X);
imx = imag(X);
modx = abs(X);
phx = atan(imx./rex);

f = fs*(0:N-1)/N;

figure;
subplot(411);
stem(f,rex);
title("rzeczywista");
subplot(412);
stem(f,imx);
title("urojona");
subplot(413);
stem(f,modx);
title("modu³");
subplot(414);
stem(f,phx);
title("faza");



















