function X =recDit1(x)

x=x(:); %w kolumnie

N = length(x);
if(N==2)    %DFT dla 2 punktów (1 motylek) i korekta e^0=1
    X(1) = x(1)+x(2);
    X(2) = x(1)-x(2);
else
   X1=recDit1(x(1:2:N));    %parzyste
   X2=recDit1(x(2:2:N));    %nieparzyste
   %X = [X1 X1] + exp(-1i*2*pi*(0:N-1)/N).*[X2 X2];  %druga polowa jest
   %symetryczna dlatego kopiuje 
   X(1:N/2) = X1 + exp(-1i*2*pi*(0:N/2-1)/N).*X2;   %'gorna czesc'
   X(N/2+1:N) = X1 + exp(-1i*2*pi*(N/2:N-1)/N).*X2;  %'dolna czesc' 
end
