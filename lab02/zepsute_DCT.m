clear all;
clc;

%--Dane--%
N = 20;
s = sqrt(1/N);

% petla do generowania wzorca cosinusowego %
for k = 1:N
    for n = 1:N
        A(k,n) = s * cos(pi*((k+0.25)-1)/N *((n-1)+0.5));
    end
    s = sqrt(2/N);
end

%Sprawdzamy czy wiersze sa ortogonalne %
nieort = 0;
for o = 1:N
    w1 = A(o,:);
    for p = 1:N
        w2 = A(p,:);
        w12 = w1 .* w2;
        if (abs(sum(w12)))>10^(-14) % ,,0'' oznacza ¿e wektory s¹ ortogonalne
            fprintf('Ta para nie jest ortogonalna: %u , %u \n', o,p);
            nieort = nieort + 1;
        end
    end
end

fprintf('Iloœæ nieortogonalnych par: %u\n', nieort);