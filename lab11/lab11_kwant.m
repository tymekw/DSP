function dq = lab11_kwant(x,Nb)
    dq = [];    
    for n = 1:2
        d = x(:,n); 
        dmax = max(d);
        dmin = min(d);
        zakres = dmax - dmin;    % minimum, maksimum, zakres
        Nq=2^Nb-1;                 % liczba przedzia³ów kwantowania
        dkw = zakres/Nq;         % szerokoœæ przedzia³u kwantowania
        dq1= dkw * round(d/dkw); % kwantyzacja sygna³u
        dq(:,n) = dq1;           % lewy i prawy kana³
    end
end