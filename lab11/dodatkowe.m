clear all; close all;
%--------------- 4 ----------------%
x1 = [0,1,2,3,3,2,1,0];
x2 = [0,7,0,2,0,2,0,7,4,2];
x3 = [0,0,0,0,0,0,0,15];

ile1 = liczbabitow(x1);
ile2 = liczbabitow(x2);
ile3 = liczbabitow(x3);
disp([num2str(ile1), ' dla x1.']);
disp([num2str(ile2), ' dla x2.']);
disp([num2str(ile3), ' dla x3.']);

%--------------- 5 ----------------%

%kodx2 = [0'110'0'10'0'10'0'110'111'10]

rng(0);
x4 = randi([1 5],1,10);
x1=x4;
    p = unique(x1);
    ile = zeros(length(p),1);
    x = sort(x1);
    for j = 1:length(p)
        for i = 1:length(x)
            if x(i) == p(j)
                ile(j) = ile(j)+1;
            end
        end
    end
    for n =1:length(ile)
        prob(n) = ile(n)/length(x1);
    end

