function H = liczbabitow(x)
    x1=x;
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
    H = -sum(prob.*log2(prob));
end