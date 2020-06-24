close all; clear all; clc;

[x, fs] = audioread('DontWorryBeHappy.wav');
N=32;
n=0:N-1;
frameLength = N;
h = sin(pi*(n+0.5)/N);
xLeft = x(:,1);
frames = floor(length(xLeft)/frameLength);
xLeft = xLeft(1:frames*frameLength);

for k=1:N/2
    for m=1:N
        A(k,m) = sqrt(4/N)*cos((2*pi/N)*(k-1+0.5)*(m-1+0.5+N/4));
    end
end
S = A';

%for i=1:frameLength:frames-1
step = frameLength/2;
Q = 200;
out = zeros(length(xLeft),1);
for i=0:step:length(xLeft)-frameLength
    sample = (h').*xLeft((i+1):i+frameLength);
    analyse = A*sample;
    quant = round(analyse*Q);
    
    synthesis = S*quant;
    synthesis = (h').*synthesis;
    out(((i+1):i+frameLength)) = out(((i+1):i+frameLength)) + synthesis;
    
    err = abs(sample-synthesis);
    %plot(sample); hold on; plot(synthesis);
     
end
out = out/Q;
figure
plot(xLeft);
figure
plot(out);
play = [out out];

error = abs(xLeft - out);
figure
plot(error);
% soundsc(play,fs);
% pause
% soundsc(x,fs);




