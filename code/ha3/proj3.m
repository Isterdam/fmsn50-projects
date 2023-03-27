clear;
close all;
load("coal_mine_disasters.mat");

%% 1c)
t_start = 1658;
t_end = 1980;

N = 25000;
burn_in = 10000;
psi = 20;
tau = T;

for d = 2:6
    theta = gamrnd(2, 1/psi); % hyperpriors
    lambda = gamrnd(2, 1/theta, 1, d);
    rho = 0.01 * ones(d, 1);

    step = (t_end-t_start)/d;
    t = t_start:step:t_end; % initial equidistant breakpoints

    breakpoints = zeros(N, length(t));

    for j = 1:burn_in % burn-in
        theta = gamrnd(2 * length(lambda) + 2, 1./(psi + sum(lambda))); % draw from p
        lambda = lambda_posterior(theta, t, tau);
    
        [~, t] = MH(lambda, t, tau, rho(1));
    end

    for j = 1:N
        theta = gamrnd(2 * length(lambda) + 2, 1./(psi + sum(lambda)));
        lambda = lambda_posterior(theta, t, tau);
    
        [acc, t] = MH(lambda, t, tau, rho(1));
        breakpoints(j, :) = t;
    end
    
    figure;
    plot(breakpoints);
    title("chain for " + (d-1) + " breakpoint(s)");
    xlabel("samples");
    ylabel("years");
end

figure;
plot(T, 1:length(T));
hold on;
xline(mean(breakpoints(:, 2)));
xline(mean(breakpoints(:, 3)));
xline(mean(breakpoints(:, 4)));
xline(mean(breakpoints(:, 5)));
xline(mean(breakpoints(:, 6)));
xlabel("years");
ylabel("accidents");
hold off;

%% 1d)

d = 6;

step = (t_end-t_start)/d;
t = t_start:step:t_end; % initial equidistant breakpoints
breakpoints = zeros(N, length(t));

rho = 0.01 * ones(d, 1);

theta_mean = zeros(psi, 1);
theta_var = zeros(psi, 1);
lambda_mean = zeros(psi, d);
lambda_var = zeros(psi, d);
t_arr = zeros(N, length(t));
t_mean = zeros(psi, length(t));
t_var = zeros(psi, length(t));

for psi = 1:50
    theta = gamrnd(2, 1/psi); % hyperpriors
    lambda = gamrnd(2, 1/theta, 1, d);
    
    theta_temp = zeros(N, 1);
    lambda_temp = zeros(N, d);

    for j = 1:burn_in % burn-in
        theta = gamrnd(2 * length(lambda) + 2, 1./(psi + sum(lambda))); % draw from p
        lambda = lambda_posterior(theta, t, tau);
    
        [~, t] = MH(lambda, t, tau, rho(1));
    end

    for j = 1:N
        theta = gamrnd(2 * length(lambda) + 2, 1./(psi + sum(lambda)));
        lambda = lambda_posterior(theta, t, tau);
    
        [acc, t] = MH(lambda, t, tau, rho(1));

        theta_temp(j) = theta;
        lambda_temp(j, :) = lambda';
        t_arr(j, :) = t;
    end
    theta_mean(psi) = mean(theta_temp);
    theta_var(psi) = var(theta_temp);
    lambda_mean(psi, :) = mean(lambda_temp);
    lambda_var(psi, :) = var(lambda_temp);
    t_mean(psi, :) = mean(t_arr);
    t_var(psi, :) = var(t_arr);
end

figure();
plot(theta_mean);
title("Theta's Dependency on psi (mean)");
xlabel("psi");
ylabel("mean");
grid on;

figure();
plot(theta_var);
title("Theta's Dependency on psi (variance)");
xlabel("psi");
ylabel("variance");
grid on;

figure();
plot(lambda_mean);
title("Lambda's Dependency on psi (mean)");
xlabel("psi");
ylabel("mean");

figure();
plot(lambda_var);
title("Lambda's Dependency on psi (variance)");
xlabel("psi");
ylabel("variance");

figure();
plot(t_mean(:, 2:end - 1));
title("Breakpoints dependency on psi (mean)");
xlabel("psi");
ylabel("mean");

figure();
plot(t_var(:, 2:end - 1));
title("Breakpoints dependency of psi (variance)");
xlabel("psi");
ylabel("variance");

%% 1e, how rho affects theta, lambda)

psi = 20;
rho = zeros(d, 30);

for i = 1:d
   rho(i, :) = (1:30) * .01; 
end

theta_mean = zeros(30, 1);
lambda_mean = zeros(30, d);

for n = 1:30
    theta = gamrnd(2, 1/psi); % hyperpriors
    lambda = gamrnd(2, 1/theta, 1, d);
    
    theta_temp = zeros(N, 1);
    lambda_temp = zeros(N, d);

    for j = 1:burn_in % burn-in
        theta = gamrnd(2 * length(lambda) + 2, 1./(psi + sum(lambda))); % draw from p
        lambda = lambda_posterior(theta, t, tau);
    
        [~, t] = MH(lambda, t, tau, rho(1, n));
    end

    for j = 1:N
        theta = gamrnd(2 * length(lambda) + 2, 1./(psi + sum(lambda)));
        lambda = lambda_posterior(theta, t, tau);
    
        [acc, t] = MH(lambda, t, tau, rho(1, n));

        theta_temp(j) = theta;
        lambda_temp(j, :) = lambda';
    end
    theta_mean(n) = mean(theta_temp);
    lambda_mean(n, :) = mean(lambda_temp);
end

figure();
plot(rho(1, :), theta_mean);
title("Theta's Dependency on rho (mean)");
xlabel("rho");
ylabel("mean");

figure();
plot(rho(1, :), lambda_mean);
title("Lambda's Dependency on rho (mean)");
xlabel("rho");
ylabel("mean");

%% 1e cont., how rho affects t)

rhos = [.01, .02, .03, .04];
rho = zeros(d, length(rhos));

for i = 1:d
   rho(i, :) = rhos; 
end

t_arr = zeros(N, length(t));

for n = 1:length(rhos)
    theta = gamrnd(2, 1/psi); % hyperpriors
    lambda = gamrnd(2, 1/theta, 1, d);

     t_arr = zeros(N, length(t));

    for j = 1:burn_in % burn-in
        theta = gamrnd(2 * length(lambda) + 2, 1./(psi + sum(lambda))); % draw from p
        lambda = lambda_posterior(theta, t, tau);
    
        [~, t] = MH(lambda, t, tau, rho(1, n));
    end

    for j = 1:N
        theta = gamrnd(2 * length(lambda) + 2, 1./(psi + sum(lambda)));
        lambda = lambda_posterior(theta, t, tau);
    
        [acc, t] = MH(lambda, t, tau, rho(1, n));

        t_arr(j, :) = t;
    end

    figure();
    subplot(2, 3, 1);
    autocorr(t_arr(:, 2), 5e2)
    title("t_1 corr func");
    xlabel("Time lag");
    ylabel("Dependency on rho = " + rho(1, n));
    
    subplot(2, 3, 2)
    autocorr(t_arr(:, 3), 5e2)
    title("t_2 corr func");
    xlabel("Time lag");
    ylabel("Dependency on rho = " + rho(1, n));
    
    subplot(2, 3, 3)
    autocorr(t_arr(:, 4), 5e2)
    title("t_3 corr func");
    xlabel("Time lag");
    ylabel("Dependency on rho = " + rho(1, n));
    
    subplot(2, 3, 4)
    autocorr(t_arr(:, 5), 5e2)
    title("t_4 corr func");
    xlabel("Time lag");
    ylabel("Dependency on rho = " + rho(1, n));

    subplot(2, 3, 5)
    autocorr(t_arr(:, 6), 5e2)
    title("t_5 corr func");
    xlabel("Time lag");
    ylabel("Dependency on rho = " + rho(1, n));

    figure();
    plot(t_arr);
    title("The breakpoints dependecy on rho, rho = " + rho(1, n));
    xlabel("Samples");
    ylabel("Year");
end

%% 1e cont., acceptance rate)

rhos = 0.001:0.001:0.1;
rho = zeros(d, length(rhos));

for i = 1:d
   rho(i, :) = rhos; 
end

acc_tot = zeros(length(rhos), d - 1);

for n = 1:length(rhos)
    theta = gamrnd(2, 1/psi); % hyperpriors
    lambda = gamrnd(2, 1/theta, 1, d);

    for j = 1:N
        theta = gamrnd(2 * length(lambda) + 2, 1./(psi + sum(lambda)));
        lambda = lambda_posterior(theta, t, tau);
    
        [acc, t] = MH(lambda, t, tau, rho(1, n));

        acc_tot(n, :) = acc_tot(n, :) + acc;
    end
end

acc_ratio = sum(acc_tot, 2) / (N * (d - 1));

figure();
plot(rhos, acc_ratio);
hold on;
plot([0, .1],[.3, .3]);
title("Acceptance ratio in different rhos");
xlabel("rho");
ylabel("acceptance ratio");
hold off;

%% 2b)

atl_txt = fopen('atlantic.txt','r');
atl_data = fscanf(atl_txt, '%f');
fclose(atl_txt);

F_inverse = @(u, beta, mu) mu - beta * log(-log(u));
n = length(atl_data);
B = 200;
[beta_hat, mu_hat] = est_gumbel(atl_data);
beta_boot = zeros(1,B);
mu_boot = zeros(1,B);
for b = 1:B % bootstrap
    y_boot = F_inverse(rand(n, 1), beta_hat, mu_hat);
    [beta, mu] = est_gumbel(y_boot);
    beta_boot(b) = beta;
    mu_boot(b) = mu;
end
beta_delta = sort(beta_boot - beta_hat); % sorting to obtain quantiles
mu_delta = sort(mu_boot - mu_hat);
alpha = 0.05; % CB level
beta_L = beta_hat - beta_delta(ceil((1 - alpha/2)*B)); % constructing CB
beta_U = beta_hat - beta_delta(ceil(alpha*B/2));
mu_L = mu_hat - mu_delta(ceil((1 - alpha/2)*B));
mu_U = mu_hat - mu_delta(ceil(alpha*B/2));
disp("95% confidence interval for beta, LB: " + beta_L + ", UB: " + beta_U);
disp("95% confidence interval for mu, LB: " + mu_L + ", UB: " + mu_U)

%% 2c)

T = 3 * 14 * 100;
wave_boot = zeros(1,B);
for b = 1:B % bootstrap
    wave_boot(b) = F_inverse(1 - 1/T, beta_boot(b), mu_boot(b));
end

wave_avg = F_inverse(1 - 1 / T, beta_hat, mu_hat);
wave_delta = sort(wave_boot(1, :) - wave_avg);
wave_U = wave_avg - wave_delta(ceil(alpha*B));
disp("95% one sided confidence interval upper bound for 100-year Atlantic wave: " + wave_U);
