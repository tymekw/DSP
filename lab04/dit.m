function X = dit( x )
N = length(x);

X1 = fft( x(1:2:N-1) ); % próbki parzyste
X2 = fft( x(2:2:N) ); % próbki nieparzyste
X = zeros( size(x) );

%k = (0:N/2-1)';
%X(1:N/2) = X1 + exp( -j*2*pi/N .*k ) .* X2;
%X(N/2+1:N) = X1 + exp( -j*2*pi/N .* (k+N/2) ) .* X2;

%to samo co w kodzie z instrukcji ale po³¹czone w 1 
X = [X1 X1] + exp(-1i*2*pi*(0:N-1)/N) .* [X2 X2];