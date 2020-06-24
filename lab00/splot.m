function y=splot(x,h)

    dlugosc=length(x)+length(h) -1;
    
    for n=1:dlugosc
        suma = 0;
        for k=1:length(h)
            if (n-k) >= 1 && n-k <= length(x)
                suma = suma + h(k)*x(n-k);  %n-k+1
            end
        end
    y(n) = suma;
end