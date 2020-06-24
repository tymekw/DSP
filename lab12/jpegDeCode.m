%function obraz = jpegDeCode( bits)
q = bi2de(bits(1:16)','left-msb');% kwant
x = bi2de(bits(17:32)','left-msb');% szerokoœæ
y = bi2de(bits(33:48)','left-msb');% wysokoœæ
bits = bits(49:end);


n=1;
prewDC = 0;
wynik = [];
while length(bits) > 4
    blok = [zeros(64,1)];
    ile_dlaDC = bi2de(bits(1:4)','left-msb'); %ile bitów dla DC
    bits = bits(5:end);
    if(ile_dlaDC == 0)
        DC = 0;
    else
        DC = bi2de(bits(1:ile_dlaDC)','left-msb'); % DC
        bits = bits(ile_dlaDC+1:end);
    end
    para1 = bi2de(bits(1:4)','left-msb'); % 1 znak z pary
    bits = bits(5:end);
    para2 = bi2de(bits(1:4)','left-msb'); % 2 znak z pary liczba bitow modulu wspolczynnika
    bits = bits(5:end);
    
    blok(1) = DC;
    if(para1 == 0 && para2 ==0)
        blok(2:end) = 0;
    end
    if(para1 ~= 0 )
        blok(2:para1) = 0;
    end
    if(para2 ~= 0)
        para3  = bi2de(bits(1:para2)','left-msb');
        if(bits(1) == 0)
            para3 = -para3;
        end
        blok(para1+1) = para3;
        bits = bits(para2+1:end);
    end
    wynik = [wynik;blok];
end
for n=1:64:length(wynik)
    blok
end


%{

bits = [];
[ xi, yi ] = size( img );
DC = 0;
for x = 1:8:xi
    for y = 1:8:yi
        blok = img( x:x+7, y:y+7 );
        dblok = DCT8x8( blok );         % realizacja dyskretnej transformaty kosinusowej
        qblok = kwant( dblok, q );      % kwantyzacja wspolczynnikow transformaty DCT
        zblok = ZigZag( qblok );        % rozwijanie bloku 8x8 do wektora 64x1 lgorytmem Zig-Zag
        [DCnew, pair] = RLE( zblok );   % tworzenie ,,par'' lgorytmem RLE
        bDC = VLCDC( DC-DCnew );        % kodowanie DC do bitow algorytmem VLC
        DC=DCnew;
        bblok = VLC( pair );            % kodowanie par do bitow algorytmem VLC
        bits = [ bits; bDC; bblok ];    % formowanie strumienia bitowego
    end
end
bq = double( dec2bin( q, 16 ) ) - 48;   % zakoduj wspolczynnik kompresji na16 bitach
[x,y] = size( img );                    % oblicz rozmiar obrazka
bx = double( dec2bin( floor(x/8)*8, 16 ) ) - 48; % zakoduj szerokosc obrazka na 16 bitach
by = double( dec2bin( floor(y/8)*8, 16 ) ) - 48; % zakoduj wysokosc obrazka na 16 bitach
bits = [ bq'; bx'; by'; bits ];         % dolaczenie na poczatku strumienia dodatkowe dane niezbedne dla dekodera


function y = DCT8x8( x )
% x: blok 8x8 pix
% y: blok 8x8 pix po transformacie DCT
  y = dct2(x);


function y = kwant( x, q )
% x: blok 8x8 wspolczynnikow transformaty DCT
% q: wspolczynnik kwantyzacji
% y: skwantowane wspolczynniki x
  y = floor(x/q + 0.5);

function y = ZigZag( x )
% x: blok 8x8 skwantowanych wspolczynnikow DCT
% y: przeksztalcenie x do wektora 64x1 algorytem ZigZag     


macierz = [0, 1, 5, 6, 14,16,27,28;
           2, 4, 7, 13,16,26,29,42;
           3, 8 ,12,17,25,30,41,43;
           9, 11,18,24,31,40,44,53;
           10,19,23,32,39,45,52,54;
           20,22,33,38,46,51,55,60;
           21,34,37,47,50,56,59,61;
           35,36,48,49,57,58,62,63];
       
[mac,idx] = sort(macierz(:));
y = [x(idx)];


function [ DC, y ] = RLE( x )
% x: skwantowane wspolczynniki w formacie wektora o wymiarach 64x1
% DC: wyliczona wartosc DC z biezacej ramki
% y: wyznaczone ,,pary'' w formacie Nx2 gdzie N oznacza ilosc par (kazdy wiersz to jedna para)
%      ilosc par w konkretnym bloku 8x8 zalezy od danych danych wejsciowych

y = [];
DC = x(1);

koniec = 64;
while x( koniec ) == 0
    koniec = koniec - 1;
end

i=2;
while i<=koniec
    if x(i)~=0
        ilezer=0;
        liczba=x(i);
        i=i+1;
    else
        ilezer=0;
        while( x(i)==0 && ilezer<15 )
            i=i+1;
            ilezer=ilezer+1;
        end
        liczba=x(i);
        i=i+1;
    end
    y =[ y; [ilezer,liczba] ];
end

if( x(64) == 0 )
    y = [ y ;[0,0] ];
end


function y = VLCDC( x )
% x: wspolczynnik DC
% y: strumien bitowy wyznaczony dla wartosci DC
%     strumien bitowy ma byc w formacie wektora Nx1 gdzie N oznacza liczbe
%     bitow, kazdy element wektora moze przyjmowac wartosc ,,0'' lub ,,1'',
%     wektor moze byc typu ,,double''

if( x == 0 )
    y = [ 0, 0, 0, 0 ]';
else
    b = double( dec2bin( abs(x) ) ) - 48;
    lb = length( b );
    if x<0
        b = 1 - b;
    end
    blb = double( dec2bin( lb,4) ) - 48;
    y = [ blb, b ]';
end

function y = VLC( x )
% x: wektor ,,par'' o wymiarach Nx2 (kazdy wiersz to jedna para)
% y: strumien bitowy wyznaczony dla wszystkich ,,par'',
%     strumien bitowy ma byc w formacie wektora Nx1 gdzie N oznacza ilosc
%     bitow, kazdy element wektora moze przyjmowac wartosc ,,0'' lub ,,1'',
%     wektor moze byc typu ,,double''

y = [];
[ xx, ~ ] = size( x );

for i=1:xx
    pair = x( i, : );
    
    xxxx = double( dec2bin( pair(1), 4 ) ) - 48;
    
    b = double( dec2bin( abs(pair(2)) ) ) - 48;
    lb = length( b );
    if pair(2)<0
        b = 1 - b;
    end
    blb = double( dec2bin( lb,4 ) ) - 48;
    if( pair(2) == 0 )
        tmp = [ xxxx, [ 0, 0, 0, 0 ] ]';
    else
        tmp = [ xxxx, blb, b ]';
    end
    y = [y; tmp];
end
%}