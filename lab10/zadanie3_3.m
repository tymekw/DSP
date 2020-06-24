% ----------------------------------------------------------
% Tabela 19-4 (str. 567)
% ï¿½wiczenie: Kompresja sygnaï¿½u mowy wedï¿½ug standardu LPC-10
% ----------------------------------------------------------

clear all; clf;
close all;

[x,fpr] = audioread('mowa1.wav');	      % wczytaj sygna³ mowy
%x = x(1:480); %A8
%t = (1:length(x))*1/fpr;
%plot(t,x); title('sygna³ mowy');            % poka¿ go

%soundsc(x,fpr);
%soundsc(x(39130:39700),fpr);	
%soundsc(x(39600:40400),fpr);


dzwieczna = 38500:39000;        % y
bezdzwieczna = 39500:40000;     % ch
between = 39000:39500;          % miedzy
%x = x(dzwieczna);                 % wybrana g³oska
wczytany_x = x;                 % nie bêdzie zmieniany


X=fft(x);                               % Fourier
f = (0:(length(X)/2))*fpr/length(X);    % oœ x
t=(1:length(x))*1/fpr;                  % oœ x
X = abs(X/length(x));                   % abs z niego
X = 2* X(1:length(X)/2+1);              % wziêcie 1 czêœci i 2*mocniej
figure;
subplot(411);plot(t,x); title('sygna³ mowy [t] //przed preemfaz¹//');            % poka¿ go
subplot(412);plot(f,X); title('widmo gestoœci mocy [Hz] //przed preemfaz¹//');



N=length(x);	  % d³ugoœæ sygna³u
Mlen=256; %240/256 normalnie/dodatkowe		  % d³ugoœæ okna Hamminga (d³ugoœæ analizowanego bloku próbek)
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

if(1)   % wykonuje lub nie preemfaze
    x = x_pre;
end

if(1)   % wykonuje okno hamminga lub nie
    x = x.*hamming(240)';
end

po_progowaniu = [];
Ps =[];
autokorelacja = [];
resztkowy = 0;
aa=[];
rresztkowy=[];
wwzm=[];
s1=[];
aproksymacja = 0;
%G£ÓWNA PÊTLA
for  nr = 1 : Nramek
    
    % pobierz kolejny fragment sygna³u
    n = 1+(nr-1)*Mstep : Mlen + (nr-1)*Mstep;
    bx = x(n);
    
    
    %low pass filter 
    lp = fir1(50,(3500/fpr),'low')
   % bx = filter(lp,1,bx);
    
    
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
    
    
    
    rr(1:Np,1)=(r(2:Np+1))';                       % wektor autokorelacji
    for m=1:Np
        R(m,1:Np)=[r(m:-1:2) r(1:Np-(m-1))];	   % zbuduj macierz autokorelacji
    end
    a=-inv(R)*rr;								   % oblicz wspó³czynniki filtra predykcji
    wzm=r(1)+r(2:Np+1)*a;						   % oblicz wzmocnienie
    H=freqz(1,[1;a]);					   		   % oblicz jego odp. czêstotliwoœciow¹
    
    offset=20; rmax=max( r(offset : Mlen) );	   % znajdŸ maksimum funkcji autokorelacji
    imax=find(r==rmax);						       % znajdŸ indeks tego maksimum
    if ( rmax > 0.35*r(1) )                        % g³oska dŸwiêczna (okresowa)/bezdŸwiêczna?
        T=imax;
        resztkowy = filter([1;a],1,bx)/wzm;
        
        %DODATKOWE
        Resztkowy = abs(fft(resztkowy)/length(resztkowy)); %fft
        
        fi = fir1(10, (900)/fpr,'low');
        
       % figure;
       % plot(abs(freqz(fi,1,fpr)))
       % pause;
       
        wygladzone = filter(fi,1,Resztkowy(1:128));
        
        
        razem = [wygladzone'; Resztkowy(129:end)'];
        %razem = Resztkowy';
        figure;
        plot(razem);
        hold on;
        plot(Resztkowy);
        pause;
        aproksymacja = polyfit(linspace(0,1,256),razem',10);
        %figure;
        %fpofil = (0:length(razem)-1)*fpr/length(razem);
        %plot(fpofil,razem','-');
        %hold on;
        %p = polyval(aproksymacja,linspace(0,1,256));
        %plot(p,'r*');
       % pause;
    else
        T=0;
    end
    % if (T>80) T=round(T/2); end				   % znaleziono drug¹ podharmoniczn¹
    T                                              % wyœwietl wartoœæ T
    if(T~=0)
        disp(["okres: ", T*1/fpr, "[s]; czestotliwoœæ: ", 1/(T*1/fpr),"[Hz]"]);
    end
    
    
    if(0)
        figure;                                        % odpowiedŸ impulsowa filra na ramkê
        f=(0:length(H)-1)*(fpr/2)/length(H);
        plot(f,abs(H)); 
        title('widmo filtra traktu g³osowego [Hz]');
    else
        lpc=[lpc; T; wzm; a; aproksymacja'];					   % KOMPRESJA!!! zapamiêtaj watoœci paramterów
    end
    
    % SYNTEZA - odtwórz na podstawie parametrów ----------------------------------------------------------------------
    
    if (T~=0)gdzie=gdzie-Mstep; end			   % "przenieœ" pobudzenie dŸwiêczne
    for n=1:Mstep
        wektor = polyval(aproksymacja,linspace(0,1,256));% wracam do wektora
        resztkowy_ap =real(ifft(wektor)*length(wektor));
        ss1 = filter(1,[1;a],resztkowy_ap)*wzm;
        if( T==0)
            pob=2*(rand(1,1)-0.5); gdzie=(3/2)*Mstep+1; % pobudzenie szumowe
        else
            if(1)%DODATKOWE
                wektor = polyval(aproksymacja,linspace(0,1,256));% wracam do wektora
                resztkowy_ap =real(ifft(wektor)*length(wektor));
                pob = resztkowy_ap(n);
                gdzie = gdzie+T;
            else
                if (n==gdzie) pob=resztkowy(n); gdzie=gdzie+T;	        % pobudzenie dŸwiêczne
                else pob=resztkowy(n); end
            end
        end
        ss(n)=wzm*pob-bs*a;		% filtracja "syntetycznego" pobudzenia
        bs=[ss(n) bs(1:Np-1) ];	% przesuniêcie bufora wyjœciowego
    end
    s = [s ss];						% zapamiêtanie zsyntezowanego fragmentu mowy
    s1 = [s1; ss1'];
end
%s = s1;
s=filter(1,[1 -0.9735],s); % filtracja (deemfaza) - filtr odwrotny - opcjonalny

 
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
else
    figure;                             % porównanie czasowe sygna³u wczytanego i odtworzonego
    t = (1:length(wczytany_x))*1/fpr;
    subplot(211); plot(t, wczytany_x);title('wczytana mowa [t]');
    t = (1:length(s))*1/fpr;
    subplot(212); plot(t, s); title('mowa zsyntezowana [t]');

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


