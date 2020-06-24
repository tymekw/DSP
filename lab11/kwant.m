function dq = kwant(x,Nb)
    dq = [];    
    d = x;
    dmax = max(d);
    dmin = min(d);
    zakres = dmax - dmin;    % minimum, maksimum, zakres
    Nq=2^Nb-1;                 % liczba przedzia��w kwantowania
    dkw = zakres/Nq;         % szeroko�� przedzia�u kwantowania
    dq= dkw * round(d/dkw);  % kwantyzacja sygna�u
end