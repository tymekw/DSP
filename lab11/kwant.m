function dq = kwant(x,Nb)
    dq = [];    
    d = x;
    dmax = max(d);
    dmin = min(d);
    zakres = dmax - dmin;    % minimum, maksimum, zakres
    Nq=2^Nb-1;                 % liczba przedziałów kwantowania
    dkw = zakres/Nq;         % szerokość przedziału kwantowania
    dq= dkw * round(d/dkw);  % kwantyzacja sygnału
end