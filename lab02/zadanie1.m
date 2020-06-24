clear all; close all;

%wzorce kosinusowe DCT-II
%n -ta próbka k-tego sygna³u (funkcji bazowej)
%funcje (sygnaly) w wierszach
N=20;
s = sqrt(1/N); %s0
for k =1:20
    for n=1:20
        A(k,n) = s*cos(pi*(k-1)/N * ((n-1)+0.5));
    end
    s = sqrt(2/N); %sk
end

%czy s¹ ortonormalne?
%iloczyn sklarny wszykich par jest równy zero?

for i=1:N
    w1 = A(i,:);
    for j=i:N   %bo bez sensu wiersz1 z wiersz2 i osobno wiersz2 z wiersz1
        w2 = A(j,:);
        iloczyn_bezposrednio(i,j) = w1*w2';
        if(abs(iloczyn_bezposrednio(i,j))<10^(-13))     %==0 
           % fprintf("ta para jest ortogonalna: %u i %u \n", i,j);
        else
            fprintf("ta para nie jest ortogonalna: %u i %u \n", i,j);
        end
        iloczyn_i_suma(i,j) = sum(w1.*w2);
    end
end


blad = iloczyn_bezposrednio-iloczyn_i_suma;
for i=1:N
for j=1:N
    if(abs(blad(i,j))>10^(-13))
        fprintf("slabo %u , %u \n",i,j)
    end
end
end

