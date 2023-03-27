function a = t_posterior(lambda, t, tau)
    n = zeros(length(t)-1, 1);
    t_diff = zeros(length(t)-1, 1);
    
    for i = 1:(length(t)-1)
       n(i) = sum((t(i) <= tau) & (tau < t(i+1)));
    end
    
    for i = 1:(length(t)-1)
        t_diff(i) = t(i+1) - t(i);
    end

    a = exp(sum(log(lambda).*n + log(t_diff) - lambda.*t_diff));
end

