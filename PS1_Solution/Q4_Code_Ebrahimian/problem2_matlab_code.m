
clear
clc
clear global

global alpha b1 b0 r lambda mu sigma

b0 = 0.1;
b1 = 1.2;
r = 0.05;
lambda = 0.02;
mu = 0.05;
sigma = 0.02;
alpha = 0.9;
k0 = 1.5;
knew = 2.5;

beta = roots([sigma^2/2, mu - sigma^2/2, - r - lambda]);
beta2 = max(beta);
clear beta
zb = (beta2/(beta2 - 1))*((r+lambda-mu)/(knew^alpha - k0^alpha))...
    *(b0 + b1*(knew - k0));

N = 4000;
T = 4000;

znew = zeros(N, T);
k = zeros(N, T);
znew(:, 1) = log(0.01);
k(:, :) = k0;


for t = 2:T
    
    znew(:, t) = znew(:, t-1) + (mu - sigma^2/2) + sigma*normrnd(0, 1, N, 1);
    l = rand(N, 1);
    znew(:, t) = znew(: , t).*(l > lambda) + log(0.01)*(l < lambda);
        
    k(:, t) = knew*(znew(:, t) > log(zb)) + ...
     knew*(k(:, t - 1) == knew).*...
     (znew(:, t) ~=log(0.01)).*(znew(:, t) < log(zb)) + ...
        k0*(k(:, t - 1) ==...
        k0).*(znew(:, t) < log(zb)).*(znew(:, t) ~=log(0.01)) + ...
                k0*(znew(:, t) == log(0.01));
    
    disp(['period ', num2str(t), ' passed']);
end


x = linspace(-5, 15, 31);
t1 = 4000;
z1 = znew(:, t1);
k1 = k(:, t1);
histogram(z1, x, 'Normalization', 'probability');
histogram(z1(k1 == k0), x, 'Normalization', 'probability');
hold on;
histogram(z1(k1 == knew), x, 'Normalization', 'probability');








