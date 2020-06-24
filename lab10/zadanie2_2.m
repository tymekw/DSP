% ----------------------------------------------------------
% Tabela 19-4 (str. 567)
% ï¿½wiczenie: Kompresja sygnaï¿½u mowy wedï¿½ug standardu LPC-10
% ----------------------------------------------------------

clear all; clf;
close all;

[x,fpr] = audioread('mowa1.wav');	      % wczytaj sygna³ mowy
[kropa4,fpr] = audioread('coldvox.wav');
%x = x(1:480); %A8 punkt 1 tak ³atwiej ?
%t = (1:length(x))*1/fpr;
%plot(t,x); title('sygna³ mowy');            % poka¿ go

%soundsc(x,fpr);
%soundsc(x(39130:39700),fpr);	
%soundsc(x(39600:40400),fpr);


dzwieczna = 38500:39000;        % y
bezdzwieczna = 39500:40000;     % ch
between = 39000:39500;          % miedzy
%x = x(dzwieczna);              % wybrana g³oska
wczytany_x = x;                 % nie bêdzie zmieniany

x_resztka = x(dzwieczna);       % punkt pierwszy?              

X=fft(x);                               % Fourier
f = (0:(length(X)/2))*fpr/length(X);    % oœ x
t=(1:length(x))*1/fpr;                  % oœ x
X = abs(X/length(x));                   % abs z niego
X = 2* X(1:length(X)/2+1);              % wziêcie 1 czêœci i 2*mocniej
figure;
subplot(411);plot(t,x); title('sygna³ mowy [t] //przed preemfaz¹//');            % poka¿ go
subplot(412);plot(f,X); title('widmo gestoœci mocy [Hz] //przed preemfaz¹//');



N=length(x);	  % d³ugoœæ sygna³u
Mlen=240;		  % d³ugoœæ okna Hamminga (d³ugoœæ analizowanego bloku próbek)
Mstep=180;		  % przesuniêcie próbek w czasie (w próbkach)
Np=10;			  % rz¹d filtra
gdzie=Mstep+1;	  % pocz¹tkowa pozycja pierwszego pobudzenia

lpc=[];		                    % tablica na wspó³czynniki modelu sygna³u mowy
s=[];				            % ca³a mowa zsyntezowana
ss=[];						    % fragment sygna³u mowy zsyntezowany
bs=zeros(1,Np);					% bufor na fragment sygna³u mowy
Nramek=floor((N-Mlen)/Mstep+1);	% ile fragmentów (ramek) jest do przetworzenia

x_pre = filter([1 -0.9735], 1, x);	% filtracja wstêpna (preemfaza) -
% opcjonalna w sumie to chyba ma wzmocniæ wysokie czêstotliwoœci

XPre = fft(x_pre);                  % Fourier po preemfazie
XPre = abs(XPre/length(x_pre));     % abs z tego znormalizowny do d³ugoœci
XPre = 2* XPre(1:length(XPre)/2+1); % bior¹ pierwsz¹ czêœæ i 2*moc

% dokoñczenie powy¿szego plotu
subplot(413);plot(t,x_pre); title('sygna³ mowy [t] //po//');            % poka¿ go
subplot(414);plot(f,XPre); title('widmo gestoœci mocy [Hz] //po//');

if(0)   % wykonuje lub nie preemfaze
    x = x_pre;
end

if(0)   % wykonuje okno hamminga lub nie
    x = x.*hamming(240)'
end

po_progowaniu = [];
Ps =[];
autokorelacja = [];


bx1 = x_resztka(50:290);
bx1 = bx1 - mean(bx1);  % usuñ wartoœæ œredni¹
for k = 0 : Mlen-1
    r1(k+1) = sum( bx1(1 : Mlen-k) .* bx1(1+k : Mlen) ); % funkcja autokorelacji
end
rr1(1:Np,1)=(r1(2:Np+1))';                       % wektor autokorelacji
for m=1:Np
  R1(m,1:Np)=[r1(m:-1:2) r1(1:Np-(m-1))];	   % zbuduj macierz autokorelacji
end
a1=-inv(R1)*rr1;								   % oblicz wspó³czynniki filtra predykcji
wzm1=r1(1)+r1(2:Np+1)*a1;						   % oblicz wzmocnienie
resztka = filter([1;a1],1,bx1)/wzm1;


%{
%kwantyzacja i porównanie
x_min=min(resztka); x_max=max(resztka); x_zakres=x_max-x_min; % minimum, maksimum, zakres
Nb=2; Nq=2^Nb; % liczba bitów, liczba przedzia³ów kwantowania
dx=x_zakres/Nq; % szerokoœæ przedzia³u kwantowania
xq=dx*round(resztka/dx); % kwantyzacja sygna³u


figure;
hold on;
plot(resztka);
plot(xq); 
%}

okres = resztka(100:220); %108:211
okres1 = interp(okres,2);

if(0)%????????
    okres1 = resztka(1:120);
    okres2 = resztka(121:240);
    okressr = (okres1+okres2)/2;
    okres1 = interp(okressr,2);
end
figure;
plot(okres1);



%G£ÓWNA PÊTLA
for  nr = 1 : Nramek
    
    % pobierz kolejny fragment sygna³u
    n = 1+(nr-1)*Mstep : Mlen + (nr-1)*Mstep;
    bx = x(n);
    
    %PROGOWANIE Z KSI¥¯KI PO POLSKU!!!!!!!
    P = 0.3*max(bx);
    bx1 = bx;
    for prog=1:length(bx1)
        if(bx1(prog)>=P)
            bx1(prog) = bx1(prog)-P;
        elseif(bx1(prog)<=-P)
            bx1(prog)=bx1(prog)+P;
        else
            bx1(prog)=0;
        end
    end
    Ps = [Ps,P];                            % zapisanie wyników P (progów)
    po_progowaniu = [po_progowaniu; bx1];   % zapisanie sygna³u po progowaniu
    
    if(0) %czy chce zeby bylo progowanie
        bx = bx1;
    end
    
    
    % ANALIZA - wyznacz parametry modelu ---clear all

    bx = bx - mean(bx);  % usuñ wartoœæ œredni¹
    for k = 0 : Mlen-1
        r(k+1) = sum( bx(1 : Mlen-k) .* bx(1+k : Mlen) ); % funkcja autokorelacji
    end
    autokorelacja = [autokorelacja;r']; % zapisuje do wektora wyniki autokorelacji dla ka¿dej ramki
   
    % figure; plot(n,bx); title('fragment sygna³u mowy');
    % figure; plot(r); title('jego funkcja autokorelacji');
    
    offset=20; rmax=max( r(offset : Mlen) );	   % znajdŸ maksimum funkcji autokorelacji
    imax=find(r==rmax);						       % znajdŸ indeks tego maksimum
    if ( rmax > 0.35*r(1) ) T=imax; else T=0; end  % g³oska dŸwiêczna (okresowa)/bezdŸwiêczna?
    % if (T>80) T=round(T/2); end				   % znaleziono drug¹ podharmoniczn¹
    T                                              % wyœwietl wartoœæ T
    if(T~=0)
        disp(["okres: ", T*1/fpr, "[s]; czestotliwoœæ: ", 1/(T*1/fpr),"[Hz]"]);
    end
    
    rr(1:Np,1)=(r(2:Np+1))';                       % wektor autokorelacji
    for m=1:Np
        R(m,1:Np)=[r(m:-1:2) r(1:Np-(m-1))];	   % zbuduj macierz autokorelacji
    end
    a=-inv(R)*rr;								   % oblicz wspó³czynniki filtra predykcji
    wzm=r(1)+r(2:Np+1)*a;						   % oblicz wzmocnienie
    H=freqz(1,[1;a]);					   		   % oblicz jego odp. czêstotliwoœciow¹
    
    if(0)
        figure;                                        % odpowiedŸ impulsowa filra na ramkê
        f=(0:length(H)-1)*(fpr/2)/length(H);
        plot(f,abs(H)); 
        title('widmo filtra traktu g³osowego [Hz]');
    else
        lpc=[lpc; T; wzm; a; ];					   % KOMPRESJA!!! zapamiêtaj watoœci paramterów
    end
    
    % SYNTEZA - odtwórz na podstawie parametrów ----------------------------------------------------------------------
    % T = 0; % punkt 1 -> taki szept
    % T = 0; % kropa 4 v
    if (T~=0) gdzie=gdzie-Mstep; end                    % "przenieœ" pobudzenie dŸwiêczne
    for n=1:Mstep
        if( T==0)
            pob=2*(rand(1,1)-0.5); gdzie=(3/2)*Mstep+1; % pobudzenie szumowe
           % pob = kropa4(n); gdzie=(3/2)*Mstep+1;     % kropa 4 -> diaboliczny robot :D
        else
            % T = T*2;                           % kropka 2 2/czest podstawowa -> robot
            % T = 80;                            % kropka 3 robot   
            if (n==gdzie) pob=okres1(n); gdzie=gdzie+T;  %dzwieczne                         
            else pob=okres1(n); end
        end
        ss(n)=wzm*pob-bs*a;		% filtracja "syntetycznego" pobudzenia
        bs=[ss(n) bs(1:Np-1) ];	% przesuniêcie bufora wyjœciowego
    end
    s = [s ss];						% zapamiêtanie zsyntezowanego fragmentu mowy
end
% s=filter(1,[1 -0.9735],s); % filtracja (deemfaza) - filtr odwrotny - opcjonalny

if(0) %1 -> nie robie kompresji; 0-> robiê
    figure;                             % sygna³ przed i po progowaniu wraz z progami
    P1 = Ps(1)*ones(length(t));
    P2 = Ps(2)*ones(length(t));
    subplot(211);plot(t,x,t,P1,t,P2); title("przed progowaniem [t]");
    t=t(1:(240*Nramek));
    P1 = Ps(1)*ones(length(t));
    P2 = Ps(2)*ones(length(t));
    subplot(212);plot(t,po_progowaniu,t,P1,t,P2);  title("po progowaniu [t]");

    figure;                             % funkcja autokorelacji co Mlen dla kolejych ramek
    n = 1:length(autokorelacja);
    plot(n,autokorelacja,n,P1,n,P2); title("autokorelacja [n]");

    figure;                             % porównanie czasowe sygna³u wczytanego i odtworzonego
    t = (1:length(wczytany_x))*1/fpr;
    subplot(211); plot(t, wczytany_x);title('wczytana mowa [t]');
    t = (1:length(s))*1/fpr;
    subplot(212); plot(t, s); title('mowa zsyntezowana [t]');
else
    disp(["d³ugoœæ sygna³u wejœciowego: ",length(wczytany_x),"d³ugoœæ wektora do wys³ania", length(lpc)]);
end
figure;                             % porównanie czêstotliwoœciowe sygan³u wczytanego i odtworzonego
XWcz = fft(wczytany_x);
f = fpr*(0:((length(XWcz))/2))/length(XWcz);
XWcz = abs(XWcz/length(wczytany_x));
XWcz = 2* XWcz(1:length(XWcz)/2+1);
%t = (1:length(wczytany_x))*1/fpr;

subplot(211); plot(f, XWcz);title('wczytana mowa [f]');
S = fft(s);
f1 = fpr*(0:((length(S))/2))/length(S);
S = abs(S/length(s));
S = 2* S(1:length(S)/2+1);
%t = (1:length(wczytany_x))*1/fpr;
subplot(212); plot(f1, S); title('mowa zsyntezowana [f]');

%pause;
%soundsc(s, fpr)


