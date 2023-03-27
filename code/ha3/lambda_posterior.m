function a = lambda_posterior(theta, t, tau)
    n = zeros(1, length(t)-1);
    
    for i = 1:(length(t)-1)
        n(i) = sum((t(i) <= tau) & (tau < t(i+1)));
    end
    
    a = gamrnd(n' + 2, 1./(theta + (t(2:end) - t(1:end-1))'));
end

