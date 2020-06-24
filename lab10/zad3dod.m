% ----------------------------------------------------------
% Tabela 19-4 (str. 567)
% �wiczenie: Kompresja sygna�u mowy wed�ug standardu LPC-10
% ----------------------------------------------------------
%commented plots and signal cutting

clear all; clf;
close all;

[x,fpr] = audioread('mowa1.wav');	      % wczytaj sygna� mowy
t = (1:length(x))*1/fpr;
%plot(x); title('sygna� mowy');            % poka� go

% Wybierz o sta�ej ampltudzie i cz�stotliwo�ci tonu podstawowego
%moje = 7822:8096;
%x=x(moje);
wczytany_x = x;
soundsc(x,fpr);
pause;
%{
t=(1:length(x))*1/fpr;
X=fft(x);
f = (0:(length(X)/2))*fpr/length(X);
X = abs(X/length(x));
X = 2* X(1:length(X)/2+1);
figure;
subplot(411);plot(x); title('sygna� mowy [t] //przed//');            % poka� go
subplot(412);plot(f,X); title('widmo gesto�ci mocy [Hz] //przed//');
%}


N=length(x);	  % d�ugo�� sygna�u
Mlen=256;		  % d�ugo�� okna Hamminga (d�ugo�� analizowanego bloku pr�bek)
Mstep=180;		  % przesuni�cie pr�bek w czasie (w pr�bkach)
Np=6;			  % rz�d filtra
gdzie=Mstep+1;	  % pocz�tkowa pozycja pierwszego pobudzenia

lpc=[];		                    % tablica na wsp�czynniki modelu sygna�u mowy
s=[];				            % ca�a mowa zsyntezowana
ss=[];						    % fragment sygna�u mowy zsyntezowany
bs=zeros(1,Np);					% bufor na fragment sygna�u mowy
Nramek=floor((N-Mlen)/Mstep+1);	% ile fragment�w (ramek) jest do przetworzenia

x_pre = filter([1 -0.9735], 1, x);	% filtracja wst�pna (preemfaza) -
% opcjonalna w sumie to chyba ma wzmocni� wysokie cz�stotliwo�ci

if(0)   %wykonuje lub nie preemfaze
    x = x_pre;
end

if(0)   %okno hamminga
    x = x.*hamming(240)';
end

%po_progowaniu = [];
%Ps =[];
autokorelacja = [];
resztkowy=[];

%G��WNA P�TLA
for  nr = 1 : Nramek

    % pobierz kolejny fragment sygna�u
    n = 1+(nr-1)*Mstep : Mlen + (nr-1)*Mstep;
    bx = x(n); 
   
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
   % Ps = [Ps,P];
   % po_progowaniu = [po_progowaniu; bx1]; 
    if(0) %czy chce zeby bylo progowanie
        bx = bx1;
    end
    
    % ANALIZA - wyznacz parametry modelu ---------------------------------------------------
    bx = bx - mean(bx);  % usu� warto�� �redni�
    for k = 0 : Mlen-1
        r(k+1) = sum( bx(1 : Mlen-k) .* bx(1+k : Mlen) ); % funkcja autokorelacji
    end
    autokorelacja = [autokorelacja;r'];
    %figure; plot(n,bx); title('fragment sygna�u mowy');
    %figure; plot(r); title('jego funkcja autokorelacji');
    
    offset=20; rmax=max( r(offset : Mlen) );	   % znajd� maksimum funkcji autokorelacji
    imax=find(r==rmax);						       % znajd� indeks tego maksimum
    if ( rmax > 0.35*r(1) ) T=imax; else T=0; end  % g�oska d�wi�czna (okresowa)/bezd�wi�czna?
    % if (T>80) T=round(T/2); end				   % znaleziono drug� podharmoniczn�
    T                                              % wy�wietl warto�� T
    if(T~=0)
       % disp(["okres: ", T, "; czestotliwo��: ", 1/T]);
       % T = T*2; %kropka 2 obnizenie g�osu, ma�o czytelne 'metaliczne'
    end
    
    rr(1:Np,1)=(r(2:Np+1))';                       % wektor autokorelacji
    for m=1:Np
        R(m,1:Np)=[r(m:-1:2) r(1:Np-(m-1))];	   % zbuduj macierz autokorelacji
    end
    a=-inv(R)*rr;								   % oblicz wsp�czynniki filtra predykcji
    wzm=r(1)+r(2:Np+1)*a;						   % oblicz wzmocnienie
    H=freqz(1,[1;a]);					   		   % oblicz jego odp. cz�stotliwo�ciow�
    
    %resztkowy1  = filter([1;a],1,wczytany_x)/wzm;
    
    %figure; 
    %f=(0:length(H)-1)*(fpr/2)/length(H);
    %plot(f,abs(H)); 
    %title('widmo filtra traktu g�osowego [Hz]');
    
    %lpc=[lpc; T; wzm; a; ];					   % zapami�taj wato�ci paramter�w
    
    % SYNTEZA - odtw�rz na podstawie parametr�w ----------------------------------------------------------------------
    % T = 0;  % kropka 1 i 4 jak pociag na torach                                     % usu� pierwszy znak "%" i ustaw: T = 80, 50, 30, 0 (w celach testowych)
    % T = 80; %kropka 3 % robot ten sam 'ton'
    %resztkowy = resztkowy1(1:Mstep);
   
    if (T~=0) gdzie=gdzie-Mstep; end			   % "przenie�" pobudzenie d�wi�czne
    for n=1:Mstep
        % T = 70; % 0 lub > 25 - w celach testowych
        if( T==0)
            pob=2*(rand(1,1)-0.5); gdzie=(3/2)*Mstep+1; % pobudzenie szumowe
        else
            resztkowy = filter([1;a],1,bx);
            pob = resztkowy(n);
        end     
        k = bs*(a);
        ss(n)=wzm*pob - bs*(a);		% filtracja "syntetycznego" pobudzenia
        bs=[ss(n) bs(1:Np-1) ];	% przesuni�cie bufora wyj�ciowego
    end
    %subplot(414); plot(ss); title('zsyntezowany fragment sygna�u mowy'); pause
    s = [s ss];						% zapami�tanie zsyntezowanego fragmentu mowy
end
%{
figure;
P1 = Ps(1)*ones(length(t));
P2 = Ps(2)*ones(length(t));
subplot(211);plot(t,x,t,P1,t,P2); title("przed progowaniem [t]");
t=t(1:240*Nramek);
P1 = Ps(1)*ones(length(t));
P2 = Ps(2)*ones(length(t));
subplot(212);plot(t,po_progowaniu,t,P1,t,P2);  title("po progowaniu [t]");

figure;
n = 1:length(autokorelacja);
plot(n,autokorelacja,n,P1,n,P2); title("autokorelacja [n]");
%}

% s=filter(1,[1 -0.9735],s); % filtracja (deemfaza) - filtr odwrotny - opcjonalny


figure;
t = (1:length(wczytany_x))*1/fpr;
subplot(211); plot(t, wczytany_x);title('wczytana mowa [t]');
t = (1:length(s))*1/fpr;
subplot(212); plot(t, s); title('mowa zsyntezowana [t]');

figure;
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

pause;
soundsc(s, fpr)

 

