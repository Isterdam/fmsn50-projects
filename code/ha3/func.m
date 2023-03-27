function res = func(lambda, t, tau)
    d = length(t) - 1;
    t_diff = zeros(d, 1);

    n = zeros(d, 1);
    
    for i = 1:d
        bla = find(tau >= t(i) & tau < t(i + 1));
        n(i) = length(bla);
    end

    for j = 1:d
        t_diff(j) = t(j + 1) - t(j);
    end
    
    res = exp(sum(log(lambda) .* n + log(t_diff) - lambda .* t_diff));
end