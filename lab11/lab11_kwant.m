function dq = lab11_kwant(x,Nb)
    dq = [];    
    for n = 1:2
        d = x(:,n); 
        dmax = max(d);
        dmin = min(d);
        zakres = dmax - dmin;    % minimum, maksimum, zakres
        Nq=2^Nb-1;                 % liczba przedzia��w kwantowania
        dkw = zakres/Nq;         % szeroko�� przedzia�u kwantowania
        dq1= dkw * round(d/dkw); % kwantyzacja sygna�u
        dq(:,n) = dq1;           % lewy i prawy kana�
    end
end