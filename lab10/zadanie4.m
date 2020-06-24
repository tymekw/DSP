% ----------------------------------------------------------
% Tabela 19-4 (str. 567)
% �wiczenie: Kompresja sygna�u mowy wed�ug standardu LPC-10
% ----------------------------------------------------------

clear all; clf;
close all;

[x,fpr] = audioread('mowa1.wav');	      % wczytaj sygna� mowy
%x = x(1:480); %A8
%t = (1:length(x))*1/fpr;
%plot(t,x); title('sygna� mowy');            % poka� go

%soundsc(x,fpr);
%soundsc(x(39130:39700),fpr);	
%soundsc(x(39600:40400),fpr);


dzwieczna = 38500:39000;        % y
bezdzwieczna = 39500:40000;     % ch
between = 39000:39500;          % miedzy
%x = x(dzwieczna);                 % wybrana g�oska
wczytany_x = x;                 % nie b�dzie zmieniany


X=fft(x);                               % Fourier
f = (0:(length(X)/2))*fpr/length(X);    % o� x
t=(1:length(x))*1/fpr;                  % o� x
X = abs(X/length(x));                   % abs z niego
X = 2* X(1:length(X)/2+1);              % wzi�cie 1 cz�ci i 2*mocniej
figure;
subplot(411);plot(t,x); title('sygna� mowy [t] //przed preemfaz�//');            % poka� go
subplot(412);plot(f,X); title('widmo gesto�ci mocy [Hz] //przed preemfaz�//');



N=length(x);	  % d�ugo�� sygna�u
Mlen=240;		  % d�ugo�� okna Hamminga (d�ugo�� analizowanego bloku pr�bek)
Mstep=180;		  % przesuni�cie pr�bek w czasie (w pr�bkach)
Np=20;			  % rz�d filtra
gdzie=Mstep+1;	  % pocz�tkowa pozycja pierwszego pobudzenia

lpc=[];		                    % tablica na wsp�czynniki modelu sygna�u mowy
s=[];				            % ca�a mowa zsyntezowana
ss=[];						    % fragment sygna�u mowy zsyntezowany
bs=zeros(1,Np);					% bufor na fragment sygna�u mowy
Nramek=floor((N-Mlen)/Mstep+1);	% ile fragment�w (ramek) jest do przetworzenia

x_pre = filter([1 -0.9735], 1, x);	% filtracja wst�pna (preemfaza) -
% opcjonalna w sumie to chyba ma wzmocni� wysokie cz�stotliwo�ci

XPre = fft(x_pre);                  % Fourier po preemfazie
XPre = abs(XPre/length(x_pre));     % abs z tego znormalizowny do d�ugo�ci
XPre = 2* XPre(1:length(XPre)/2+1); % bior� pierwsz� cz�� i 2*moc

% doko�czenie powy�szego plotu
subplot(413);plot(t,x_pre); title('sygna� mowy [t] //po//');            % poka� go
subplot(414);plot(f,XPre); title('widmo gesto�ci mocy [Hz] //po//');

if(1)   % wykonuje lub nie preemfaze
    x = x_pre;
end

if(1)   % wykonuje okno hamminga lub nie
    x = x.*hamming(240)';
end

po_progowaniu = [];
Ps =[];
autokorelacja = [];
aa = [];
aaq=[]
%G��WNA P�TLA
for  nr = 1 : Nramek
    
    % pobierz kolejny fragment sygna�u
    n = 1+(nr-1)*Mstep : Mlen + (nr-1)*Mstep;
    bx = x(n);
    
    %low pass filter 
    lp = fir1(48,(900/fpr),'low')
    bx = filter(lp,1,bx);
    
    %PROGOWANIE Z KSI��KI PO POLSKU!!!!!!!
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
    Ps = [Ps,P];                            % zapisanie wynik�w P (prog�w)
    po_progowaniu = [po_progowaniu; bx1];   % zapisanie sygna�u po progowaniu
    
    if(1) %czy chce zeby bylo progowanie
        bx = bx1;
    end
    
    
    % ANALIZA - wyznacz parametry modelu ---clear all
    
    bx = bx - mean(bx);  % usu� warto�� �redni�
    for k = 0 : Mlen-1
        r(k+1) = sum( bx(1 : Mlen-k) .* bx(1+k : Mlen) ); % funkcja autokorelacji
    end
    autokorelacja = [autokorelacja;r']; % zapisuje do wektora wyniki autokorelacji dla ka�dej ramki
   
    % figure; plot(n,bx); title('fragment sygna�u mowy');
    % figure; plot(r); title('jego funkcja autokorelacji');
    
    offset=20; rmax=max( r(offset : Mlen) );	   % znajd� maksimum funkcji autokorelacji
    imax=find(r==rmax);						       % znajd� indeks tego maksimum
    if ( rmax > 0.35*r(1) ) T=imax; else T=0; end  % g�oska d�wi�czna (okresowa)/bezd�wi�czna?
    % if (T>80) T=round(T/2); end				   % znaleziono drug� podharmoniczn�
    T                                              % wy�wietl warto�� T
    if(T~=0)
        disp(["okres: ", T*1/fpr, "[s]; czestotliwo��: ", 1/(T*1/fpr),"[Hz]"]);
    end
    
    rr(1:Np,1)=(r(2:Np+1))';                       % wektor autokorelacji
    for m=1:Np
        R(m,1:Np)=[r(m:-1:2) r(1:Np-(m-1))];	   % zbuduj macierz autokorelacji
    end
    a=-inv(R)*rr;								   % oblicz wsp�czynniki filtra predykcji
    wzm=r(1)+r(2:Np+1)*a;						   % oblicz wzmocnienie
    H=freqz(1,[1;a]);					   		   % oblicz jego odp. cz�stotliwo�ciow�
    
    %DODATKOWE
    %[gamma,e] = levinson(bx,Np);
    %liczenie gamma z a
    P = length(a);
    ax(P,1:P) = a(1:P);
    for i =P:-1:1
        g(i)=-ax(i,i);
        for j=1:i-1
            ax(i-1,j) = (ax(i,j)+g(i)*ax(i,i-j))/(1-g(i)^2);
        end
    end
    gamma =g;
    
    %kwantyzacja
    gaMin=min(gamma);
    gaMax=max(gamma); 
    gaZakres=gaMax-gaMin; % minimum, maksimum, zakres
    Nb=32;       % liczba bit�w %8 bez LP ; %10 z ale 16 super
    Nq=2^Nb; % liczba przedzia��w kwantowania
    dga=gaZakres/Nq; % szeroko�� przedzia�u kwantowania
    gaq=dga*round(gamma/dga); % kwantyzacja sygna�u
    
    %skwantowane gamma
    g = gaq;
    
    %kwantyzacja i por�wnanie
    aMin=min(a);
    aMax=max(a); 
    aZakres=aMax-aMin; % minimum, maksimum, zakres
    Nb=32;       % liczba bit�w %8 bez LP ; %10 z ale 16 super
    Nq=2^Nb; % liczba przedzia��w kwantowania
    da=aZakres/Nq; % szeroko�� przedzia�u kwantowania
    aq=da*round(a/da); % kwantyzacja sygna�u
        
    %aa = [aa;a];
    %aaq = [aaq;aq];
    
    
    %a=aq;
    
    if(0)
        figure;                                        % odpowied� impulsowa filra na ramk�
        f=(0:length(H)-1)*(fpr/2)/length(H);
        plot(f,abs(H)); 
        title('widmo filtra traktu g�osowego [Hz]');
    else
        lpc=[lpc; T; wzm; a];					   % KOMPRESJA!!! zapami�taj wato�ci paramter�w
    end
    
    % SYNTEZA - odtw�rz na podstawie parametr�w ----------------------------------------------------------------------
    %a=aq;
    if(0)
    if (T~=0) gdzie=gdzie-Mstep; end			   % "przenie�" pobudzenie d�wi�czne
    for n=1:Mstep
        if( T==0)
            pob=2*(rand(1,1)-0.5); gdzie=(3/2)*Mstep+1; % pobudzenie szumowe
        else
            if (n==gdzie) pob=1; gdzie=gdzie+T;	        % pobudzenie d�wi�czne
            else pob=0; end
        end
        ss(n)=wzm*pob-bs*aq;		% filtracja "syntetycznego" pobudzenia
        bs=[ss(n) bs(1:Np-1) ];	% przesuni�cie bufora wyj�ciowego
    end
    s = [s ss];						% zapami�tanie zsyntezowanego fragmentu mowy
    else
        e1 = zeros(1,P);
        for n=1:Mstep
            e0(1) = bx(n);
            for k=1:P
                e0(k+1) = e0(k)*g(k)*e1(k);
                e1n(k+1) = -g(k)*e0(k)+e1(k);
            end
            e1 = [bx(n) e1n(2:P)];
            y1(n) = e0(P+1);
        end
        s = [s y1];	
    end
    
end
% s=filter(1,[1 -0.9735],s); % filtracja (deemfaza) - filtr odwrotny - opcjonalny

if(0) %1 -> nie robie kompresji; 0-> robi�
    figure;                             % sygna� przed i po progowaniu wraz z progami
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
else
    figure;                             % por�wnanie czasowe sygna�u wczytanego i odtworzonego
    t = (1:length(wczytany_x))*1/fpr;
    subplot(211); plot(t, wczytany_x);title('wczytana mowa [t]');
    t = (1:length(s))*1/fpr;
    subplot(212); plot(t, s); title('mowa zsyntezowana [t]');

    disp(["d�ugo�� sygna�u wej�ciowego: ",length(wczytany_x),"d�ugo�� wektora do wys�ania", length(lpc)]);
end
figure;                             % por�wnanie cz�stotliwo�ciowe sygan�u wczytanego i odtworzonego
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

figure;
plot(aa);
hold on;
plot(aaq);

%pause;
soundsc(s, fpr)


