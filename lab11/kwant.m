function dq = kwant(x,Nb)
    dq = [];    
    d = x;
    dmax = max(d);
    dmin = min(d);
    zakres = dmax - dmin;    % minimum, maksimum, zakres
    Nq=2^Nb-1;                 % liczba przedzia³ów kwantowania
    dkw = zakres/Nq;         % szerokoœæ przedzia³u kwantowania
    dq= dkw * round(d/dkw);  % kwantyzacja sygna³u
end