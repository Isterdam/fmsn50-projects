function [acc, t] = MH(lambda, t, tau, rho)
    acc = zeros(1, length(t)-2);
    
    for i=2:(length(t)-1)
       R = rho*(t(i+1)-t(i-1));
       eps = unifrnd(-R, R);
       t_star = t(i) + eps;
       
       while (t_star < t(i-1) || t_star > t(i+1))
           eps = unifrnd(-R, R);
           t_star = t(i) + eps;
       end
       
       num = t_posterior(lambda, [t(1:i-1), t_star, t(i+1:end)], tau);
       den = t_posterior(lambda, t, tau);
       alpha = min(1, num/den);
       
       if (rand <= alpha)
           t(i) = t_star;
           acc(i-1) = acc(i-1)+1;
       end
    end
end

