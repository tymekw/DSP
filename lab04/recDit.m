function X =recDit(x)

N = length(x);

if(N==2)
    X(1) = x(1)+x(2);
    X(2) = x(1)-x(2);
else
   X1=recDit(x(1:2:N));
   X2=recDit(x(2:2:N));   
   X = [X1 X1] + exp(-1i*2*pi*(0:N-1)/N).*[X2 X2];
end
