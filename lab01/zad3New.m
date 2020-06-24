%Program do zadania 3 wyszukuj¹cy prefixy
clear all; close all;
%Bierzemy kolejne prefixy adsl maj¹ce 32 bity i liczymy autokorelacjê 
%Jak bêd¹ 2 piki, to znaczy, ¿e ci¹g bitów siê powtarza ergo jest prefixem
load("adsl_x.mat")

for i = 0:1985
   a = x(1+i:32+i);
   corr_vect = abs(xcorr(x, a));
   [p, j] = max(corr_vect);
   corr_vect(j ) = 0;
   [p1, j1] = max(corr_vect);
   if p1 == p
      disp(["1 bit prefixu = ", j-2048]) 
      disp(["1 bit koñca ramki = ", j1-2048]) 
   end
     %plot(abs(xcorr(x,a))); pause;
%    title(1+i)
end
load("adsl_x.mat")%Plik ³aduje siê jako x, wykorzystujemy pó¿niej tê zmienn¹
a = x(1+495:32+495); %prefix adsl
% a = x(1+1039:32+1039);
% a = x(1+1803:32+1803); %Something is no yes
%a = x(1+551:32+551);
%Wykonujemy funkcjê korelacji sygna³u, ¿eby znaleŸæ koniec ramki
corr_vect = xcorr(x, a);
%Rysujemy wykres korelacji
figure;
stem(abs(corr_vect));
%Szukamy indeksu, w którym korelacja ma peak 
[p, i] = max(corr_vect);