% 1 N-punktowa transformacja R sygna³u 
% za pomoc¹  N/2 punktowej transformaty zespolonej
clc;
clear all;
N=1024;
x = randn(1,N);

X0 = fft(x);

%krok 0
for n=1:(N/2)
   y(n) = x(2*n-1) + 1i*x(2*n);
end

Y = fft(y);

%krok 1
for k=1:N/2-1
    X(k) = 0.5*(Y(k+1)+conj(Y(N/2-k+1)))+0.5*1i*exp(-1i*2*pi*k/N)*(conj(Y(N/2-k+1))-Y(k+1));
end

X = [0, X];
X(1) = real(Y(1))+imag(Y(1));


X2 = conj(X(2:end));
for i = 1:length(X2)/2
    tmp = X2(i);
    X2(i) = X2(length(X2)-(i-1));
    X2(length(X2)-(i-1)) = tmp;
end
X2 = [0, X2];
X2(1) = real(Y(1))-imag(Y(1));



X3 = [X,X2];
figure;
stem(abs(X3),'b');
xlim([1 N]);
hold on;
stem(abs(X0));

mean(abs(X0-X3))    %blad













