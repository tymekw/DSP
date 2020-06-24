function Y = fftshift2D( X )
% przestawianie æwiartek widma 2D DFT
[M N] = size(X);
Y(M/2+1:M,N/2+1:N) = X(1:M/2,1:N/2);
Y(1:M/2,1:N/2) = X(M/2+1:M,N/2+1:N);
Y(M/2+1:M,1:N/2) = X(1:M/2,N/2+1:N);
Y(1:M/2,N/2+1:N) = X(M/2+1:M,1:N/2);