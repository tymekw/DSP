clc;
clear all;

x = [zeros(1,9) ones(1,6) zeros(1,9)]; %impuls prostokatny
h = x; %sygnal filtrujacy

y = splot(x,h);
a=1;
b=h;

y2 = filter(b,a,x);

y1 = conv(x,h);
%parameters
%'full' -> default
%'same' -> same length as x
%'valid' -> only those parts that are without zero -padded edges
% returns 6



fileID = fopen('x_h_sign.txt','w');
A = [x;h];
fprintf(fileID, '%d %d\n',A);
fclose(fileID);

[wynik] = textread('x_h_conv.txt', '%d');
wynik = wynik';



plot(y,'blue');
hold on;
plot(y1, 'red');
plot(y2, 'black');
plot(wynik,'green*');
legend('my function','conv', 'filter', 'cpp');