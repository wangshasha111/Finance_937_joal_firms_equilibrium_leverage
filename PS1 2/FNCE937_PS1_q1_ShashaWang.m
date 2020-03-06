% Shasha Wang PS1 Q1 FNCE 937
clear all;
close all;

cd 'E:\Dropbox\fall 19-20\Finance 937\PS1\question_1'
tic;

%% Q1 Solving the Problem of the Firm in Discrete Time
% Grid for k - centered around k_mean = 1. 
% k_max = (a_max/delta)^(1/(1-alpha))
% k_min = 0
% Grid for a - choose a_mean s.t. k_mean = 1. Nine grid points work well.
% k_mean = (alpha * a_mean / (r_mean + delta))^(1/(1-alpha))
% a_min = a_mean - m*sigma/(1-rho^2)^0.5
% a_max = a_mean + m*sigma/(1-rho^2)^0.5

%% Parametization
rho = 0.8;
r = 0.02;
delta = 0.1;
theta_1 = 0.3;
theta_2 = 0.6;
W = 2;
sigma = 0.1;
m = 4; % number of a is 2 * 4 + 1 = 9
Na = m*2 + 1;
Nk_1 = 200;
Nk_2 = 25;
b_0 = 0;
b_1 = 0.5;

%% Create grid for k and a
k_min = 0.01; % setting k_min to 0 would pop up errors in the iteration
% a_mean = (r + delta) / theta_1; % Lecture 1 Page 15/52 "picking E[a] s.t. k_mean=1"
% But that would make a too small
a_mean = 2; % I set the mean of a higher
log_a_mean = log(a_mean);
log_a_min=log(a_mean)-m*sigma/sqrt(1-rho^2);
log_a_max=log(a_mean)+m*sigma/sqrt(1-rho^2);
log_a_grid=linspace(log_a_min,log_a_max,Na);
a_min=exp(log_a_min);
a_max=exp(log_a_max);
%a_min = max(0.2, a_mean - m * sigma/(1-rho^2)^0.5);
%a_mean = a_min + m*sigma/(1-rho^2)^0.5; % update a_mean
k_mean = (theta_1 * a_mean / (r + delta))^(1/(1-theta_1)); % Lecture 1 Page 15/52
%a_max = a_mean + m*sigma/(1-rho^2)^0.5;
%a_grid = linspace(a_min, a_max, Na);
k_max = (exp(log_a_max)/delta)^(1/(1-theta_1));
k_grid_1 = linspace(k_min, k_max, Nk_1);
k_grid_2 = linspace(k_min, k_max, Nk_2);

%% a) Plot the profit function as a function of the capital stock
labor = @(a,k) (k.^theta_1 * theta_2 * a/W).^(1/(1-theta_2));
labor_min = labor(a_min,k_grid_1);
labor_mean = labor(a_mean,k_grid_1);
labor_max = labor(a_max,k_grid_1);

profit = @(labor,a,k)a * k.^theta_1.* labor.^theta_2 - W*labor;
investment = @(k,k_prime)k_prime - (1-delta)*k;
adjustment_cost = @(investment,k)b_0*k + b_1*k*(investment/k-delta).^2;
divident = @(profit,investment,adjustment_cost)profit-investment-adjustment_cost;

profit_min = profit(labor_min,a_min,k_grid_1);
profit_mean = profit(labor_mean,a_mean,k_grid_1);
profit_max = profit(labor_max,a_max,k_grid_1);

figure (1)
plot(k_grid_1, profit_min,k_grid_1, profit_mean, k_grid_1, profit_max);
title('Profit under different shocks')
xlabel('Capital Stock')
ylabel('Profit')
legend('a lowest', 'a mean', 'a highest')
xlim([0,k_max])
savefig('q1a')

%% a's transition matrix

[Z,a_prob] = tauchen(Na,log_a_mean,rho,sigma,4);
assert(sum(abs(Z-log_a_grid))<0.01);

a_grid = exp(log_a_grid);

%% Value Function Iteration and q1c Policy Function
%% See if number of capital grid is 200
Nk=Nk_1;

k_policy_index = zeros(Na,Nk);
k_policy = zeros(Na,Nk);
investment_policy = zeros(Na,Nk);
financing_policy = zeros(Na,Nk);
value   = zeros(Na,Nk);
value0 = ones(Na,Nk); % initial guess
tolerance = 0.00001;
iteration = 0;
distance = 100;

k_grid = linspace(k_min, k_max, Nk);
[aa,kk]=meshgrid(a_grid, k_grid);

figure(2)
mesh(aa, kk, value0')

tic;
% phi = b_0 * k_grid_1 + b_1 * (i/k_grid_1 - delta).^2 .* k_grid_1;
while distance > tolerance
    for i=1:Na
        for j=1:Nk
            a = a_grid(i);
            k = k_grid(j);
            l = labor(a,k);
            %l = (k^theta_1 * theta_2 * a/W).^(1/(1-theta_2));
            pi = profit(l,a,k);
            inv = investment(k,k_grid);
            %pi = a*k^theta_1*l^theta_2 - W*l;
            phi = adjustment_cost(inv,k);
            d = divident(pi,inv,phi);
            %d = pi + (1-delta-b_0)*k-k_grid-b_1*k*(k_grid/k-1).^2;
            [v,ind] = max(d + 1/(1+r)*a_prob(i,:)*value0);
            value(i,j)=v;
            k_policy_index(i,j)=ind;
            k_policy(i,j) = k_grid(k_policy_index(i,j));
        end
    end
    hold on
    mesh(aa, kk, value')
    distance = sum(sum(abs(value-value0)));
    value0=value;
    iteration = iteration + 1;
end
display("iteration =    " + iteration + "   difference =   " + distance)
toc;

title(['Value Function Iteration Nk=', num2str(Nk)])
xlabel('Productivity')
ylabel('Capital Stock')
zlabel('Value')
savefig('q1b1_ValueFunction_3d') 

figure(3)
plot(k_grid,value(1,:))
xlim([0,k_max])
iter = 1;
while iter<=Na
    hold on; plot(k_grid,value(iter,:));
    iter=iter+1;
end
title(['value function at different shocks Nk= ', num2str(Nk)])
xlabel('Capital Stock')
ylabel('Value')
savefig('q1b1_ValueFunction_2d') 

figure(4)
[aa,kk]=meshgrid(a_grid, k_grid);
mesh(aa, kk, k_policy')
title(['Policy Function for Capital Nk=', num2str(Nk)])
xlabel('Productivity')
ylabel('Capital Stock')
zlabel('Capital Next Period')
savefig('q1b1_k_policy_3d') 

figure(5)
plot(k_grid, k_policy(1,:), k_grid, k_policy(m+1,:), k_grid, k_policy(Na,:) )
legend("lowest a","average a","highest a")
title(['Policy Function for Next Period Capital Nk=', num2str(Nk)])
xlabel('Current Capital Stock')
ylabel('Capital Next Period')
xlim([0,k_max])
savefig('q1b1_k_policy_2d') 

save('nk200_results','k_grid','k_policy','value','Nk','k_policy_index')

%% See if number of capital grid is 25
close all;
Nk=Nk_2;

k_policy_index = zeros(Na,Nk);
k_policy = zeros(Na,Nk);
value   = zeros(Na,Nk);
value0 = ones(Na,Nk); % initial guess
tolerance = 0.00001;
iteration = 0;
distance = 100;

k_grid = linspace(k_min, k_max, Nk);
[aa,kk]=meshgrid(a_grid, k_grid);

figure(2)
mesh(aa, kk, value0')

tic;
% phi = b_0 * k_grid_1 + b_1 * (i/k_grid_1 - delta).^2 .* k_grid_1;
while distance > tolerance
    for i=1:Na
        for j=1:Nk
            a = a_grid(i);
            k = k_grid(j);
            l = labor(a,k);
            %l = (k^theta_1 * theta_2 * a/W).^(1/(1-theta_2));
            pi = profit(l,a,k);
            inv = investment(k,k_grid);
            %pi = a*k^theta_1*l^theta_2 - W*l;
            phi = adjustment_cost(inv,k);
            d = divident(pi,inv,phi);
            %d = pi + (1-delta-b_0)*k-k_grid-b_1*k*(k_grid/k-1).^2;
            [v,ind] = max(d + 1/(1+r)*a_prob(i,:)*value0);
            value(i,j)=v;
            k_policy_index(i,j)=ind;
            k_policy(i,j) = k_grid(k_policy_index(i,j));
        end
    end
    hold on
    mesh(aa, kk, value')
    distance = sum(sum(abs(value-value0)));
    value0=value;
    iteration = iteration + 1;
end
display("iteration =    " + iteration + "   difference =   " + distance)
toc;

title(['Value Function Iteration Nk=', num2str(Nk)])
xlabel('Productivity')
ylabel('Capital Stock')
zlabel('Value')
savefig('q1b2_ValueFunction_3d') 

figure(3)
plot(k_grid,value(1,:))
xlim([0,k_max])
iter = 1;
while iter<=Na
    hold on; plot(k_grid,value(iter,:));
    iter=iter+1;
end
title(['value function at different shocks Nk= ', num2str(Nk)])
xlabel('Capital Stock')
ylabel('Value')
savefig('q1b2_ValueFunction_2d') 

figure(4)
[aa,kk]=meshgrid(a_grid, k_grid);
mesh(aa, kk, k_policy')
title(['Policy Function for Capital Nk=', num2str(Nk)])
xlabel('Productivity')
ylabel('Capital Stock')
zlabel('Capital Next Period')
savefig('q1b2_k_policy_3d') 

figure(5)
plot(k_grid, k_policy(1,:), k_grid, k_policy(m+1,:), k_grid, k_policy(Na,:) )
legend("lowest a","average a","highest a")
title(['Policy Function for Next Period Capital Nk=', num2str(Nk)])
xlabel('Current Capital Stock')
ylabel('Capital Next Period')
xlim([0,k_max])
savefig('q1b2_k_policy_2d') 

save('nk25_results','k_grid','k_policy','value','Nk','k_policy_index')

%% Calculate and plot optimal investment i(a,k) and financing (-d(a,k) when positive) policies
% against the current stock of capital. Do this for both the highest,
% average, and lowest value of the shock a
% i(a,k)=k'(a,k)-(1-delta)*k
close all;

load('nk200_results.mat')

investment_policy=ones(Na, Nk);
for i=1:Na
    for j=1:Nk
        investment_policy(i,j) = investment(k_grid(j),k_policy(i,j));
    end
end

[aa,kk]=meshgrid(a_grid, k_grid);
figure(1)
mesh(aa', kk', investment_policy)
title(['Policy Function for Investment Nk=', num2str(Nk)])
xlabel('Productivity')
ylabel('Capital Stock')
zlabel('Investment')
savefig('q1c1_investment_policy_3d') 

figure(2)
plot(k_grid, investment_policy(1,:), k_grid, investment_policy(m+1,:), k_grid, investment_policy(Na,:) )
legend("lowest a","average a","highest a")
title(['Policy Function for Investment Nk=', num2str(Nk)])
xlabel('Current Capital Stock')
ylabel('Optimal Investment')
savefig('q1c1_investment_policy_2d') 

financing_policy=ones(Na, Nk); % -d(a,k) when positive
for i=1:Na
    for j=1:Nk
        a=a_grid(i);
        k=k_grid(j);
        k_prime=k_policy(i,j);
        l=labor(a,k);
        pi=profit(l,a,k);
        inv=investment_policy(i,j);
        %inv=investment(k,k_prime);
        %assert(investment==investment(i,j));
        phi=adjustment_cost(inv,k);
        d=divident(pi,inv,phi);
        financing_policy(i,j) = max(-d,0);
    end
end

figure(3)
[aa,kk]=meshgrid(a_grid, k_grid);
mesh(aa', kk', financing_policy)
title(['Policy Function for Financing Nk=', num2str(Nk)])
xlabel('Productivity')
ylabel('Capital Stock')
zlabel('Financing')
savefig('q1c1_financing_policy_3d') 

figure(4)
plot(k_grid, financing_policy(1,:), k_grid, financing_policy(m+1,:), k_grid, financing_policy(Na,:) )
legend("lowest a","average a","highest a")
title(['Policy Function for Financing Nk=', num2str(Nk)])
xlabel('Current Capital Stock')
ylabel('Optimal Financing - Max(-d,0)')
%xlim([0,20])
savefig('q1c1_financing_policy_2d') 

save('nk200_results','k_grid','k_policy','value','Nk','k_policy_index','investment_policy','financing_policy')

%% see if # of grid is 25
close all;

load('nk25_results.mat')

investment_policy=ones(Na, Nk);
for i=1:Na
    for j=1:Nk
        investment_policy(i,j) = investment(k_grid(j),k_policy(i,j));
    end
end

[aa,kk]=meshgrid(a_grid, k_grid);
figure(1)
mesh(aa', kk', investment_policy)
title(['Policy Function for Investment Nk=', num2str(Nk)])
xlabel('Productivity')
ylabel('Capital Stock')
zlabel('Investment')
savefig('q1c2_investment_policy_3d') 

figure(2)
plot(k_grid, investment_policy(1,:), k_grid, investment_policy(m+1,:), k_grid, investment_policy(Na,:) )
legend("lowest a","average a","highest a")
title(['Policy Function for Investment Nk=', num2str(Nk)])
xlabel('Current Capital Stock')
ylabel('Optimal Investment')
savefig('q1c2_investment_policy_2d') 

financing_policy=ones(Na, Nk); % -d(a,k) when positive
for i=1:Na
    for j=1:Nk
        a=a_grid(i);
        k=k_grid(j);
        k_prime=k_policy(i,j);
        l=labor(a,k);
        pi=profit(l,a,k);
        inv=investment_policy(i,j);
        %inv=investment(k,k_prime);
        %assert(investment==investment(i,j));
        phi=adjustment_cost(inv,k);
        d=divident(pi,inv,phi);
        financing_policy(i,j) = max(-d,0);
    end
end

figure(3)
[aa,kk]=meshgrid(a_grid, k_grid);
mesh(aa', kk', financing_policy)
title(['Policy Function for Financing Nk=', num2str(Nk)])
xlabel('Productivity')
ylabel('Capital Stock')
zlabel('Financing')
savefig('q1c2_financing_policy_3d') 

figure(4)
plot(k_grid, financing_policy(1,:), k_grid, financing_policy(m+1,:), k_grid, financing_policy(Na,:) )
legend("lowest a","average a","highest a")
title(['Policy Function for Financing Nk=', num2str(Nk)])
xlabel('Current Capital Stock')
ylabel('Optimal Financing - Max(-d,0)')
%xlim([0,20])
savefig('q1c2_financing_policy_2d') 

save('nk25_results','k_grid','k_policy','value','Nk','k_policy_index','investment_policy','financing_policy')

toc;

%% q1d
close all;

%% Nk=200

load('nk200_results.mat')
optimal_investment_at_mean_shock=ones(4,Nk);

%% q1d1 b_0=0,b_1=0.5
% This is just what we've been doing before
optimal_investment_at_mean_shock(1,:)=investment_policy(m+1,:);

%% q1d2 b_0=0,b_1=10
b_0=0;
b_1=10;
adjustment_cost = @(investment,k)b_0*k + b_1*k*(investment/k-delta).^2;

k_policy_index = zeros(Na,Nk);
k_policy = zeros(Na,Nk);
investment_policy2 = zeros(Na,Nk);
value   = zeros(Na,Nk);
value0 = ones(Na,Nk); % initial guess
tolerance = 0.00001;
iteration = 0;
distance = 100;

k_grid = linspace(k_min, k_max, Nk);
[aa,kk]=meshgrid(a_grid, k_grid);

tic;
% phi = b_0 * k_grid_1 + b_1 * (i/k_grid_1 - delta).^2 .* k_grid_1;
while distance > tolerance
    for i=1:Na
        for j=1:Nk
            a = a_grid(i);
            k = k_grid(j);
            l = labor(a,k);
            %l = (k^theta_1 * theta_2 * a/W).^(1/(1-theta_2));
            pi = profit(l,a,k);
            inv = investment(k,k_grid);
            %pi = a*k^theta_1*l^theta_2 - W*l;
            phi = adjustment_cost(inv,k);
            d = divident(pi,inv,phi);
            %d = pi + (1-delta-b_0)*k-k_grid-b_1*k*(k_grid/k-1).^2;
            [v,ind] = max(d + 1/(1+r)*a_prob(i,:)*value0);
            value(i,j)=v;
            k_policy_index(i,j)=ind;
            k_policy(i,j) = k_grid(k_policy_index(i,j));
            investment_policy2(i,j) = investment(k,k_policy(i,j));
        end
    end
    distance = sum(sum(abs(value-value0)));
    value0=value;
    iteration = iteration + 1;
end
display("iteration =    " + iteration + "   difference =   " + distance)
toc;

optimal_investment_at_mean_shock(2,:)=investment_policy2(m+1,:);

%plot(k_grid,optimal_investment_at_mean_shock(1,:),k_grid,optimal_investment_at_mean_shock(2,:))

%% q1d3 b_0=0,b_1=b+=0.5 when i>delta*k and b_1=b-=10 when i<delta*k
b_0=0;
b_1_plus=0.5;
b_1_minus=10;
adjustment_cost = @(investment,k)b_0*k + k*(investment/k-delta)^2*b_1_plus*(investment>delta*k)...
                                       + k*(investment/k-delta)^2*b_1_minus*(investment<=delta*k);

k_policy_index = zeros(Na,Nk);
k_policy = zeros(Na,Nk);
investment_policy3 = zeros(Na,Nk);
value   = zeros(Na,Nk);
value0 = ones(Na,Nk); % initial guess
tolerance = 0.00001;
iteration = 0;
distance = 100;

k_grid = linspace(k_min, k_max, Nk);
[aa,kk]=meshgrid(a_grid, k_grid);

tic;
% phi = b_0 * k_grid_1 + b_1 * (i/k_grid_1 - delta).^2 .* k_grid_1;
while distance > tolerance
    for i=1:Na
        for j=1:Nk
            a = a_grid(i);
            k = k_grid(j);
            l = labor(a,k);
            %l = (k^theta_1 * theta_2 * a/W).^(1/(1-theta_2));
            pi = profit(l,a,k);
            inv = investment(k,k_grid);
            %pi = a*k^theta_1*l^theta_2 - W*l;
            phi = ones(1,Nk);
            for ss=1:Nk
                phi(ss) = adjustment_cost(inv(ss),k);
            end
            d = divident(pi,inv,phi);
            %d = pi + (1-delta-b_0)*k-k_grid-b_1*k*(k_grid/k-1).^2;
            [v,ind] = max(d + 1/(1+r)*a_prob(i,:)*value0);
            value(i,j)=v;
            k_policy_index(i,j)=ind;
            k_policy(i,j) = k_grid(k_policy_index(i,j));
            investment_policy3(i,j) = investment(k,k_policy(i,j));
        end
    end
    distance = sum(sum(abs(value-value0)));
    value0=value;
    iteration = iteration + 1;
end
display("iteration =    " + iteration + "   difference =   " + distance)
toc;

optimal_investment_at_mean_shock(3,:)=investment_policy3(m+1,:); 

%% q1d4 b_0=0.02,b_1=0.5 unless i=delta*k, i.e., just replacing...
% depreciated capital takes no adjustment cost
b_0=0.02;
b_1=0.5;
adjustment_cost = @(investment,k)b_0*k*(investment~=delta*k) + b_1*k*(investment/k-delta).^2;

k_policy_index = zeros(Na,Nk);
k_policy = zeros(Na,Nk);
investment_policy4 = zeros(Na,Nk);
value   = zeros(Na,Nk);
value0 = ones(Na,Nk); % initial guess
tolerance = 0.00001;
iteration = 0;
distance = 100;

k_grid = linspace(k_min, k_max, Nk);
[aa,kk]=meshgrid(a_grid, k_grid);

tic;
% phi = b_0 * k_grid_1 + b_1 * (i/k_grid_1 - delta).^2 .* k_grid_1;
while distance > tolerance
    for i=1:Na
        for j=1:Nk
            a = a_grid(i);
            k = k_grid(j);
            l = labor(a,k);
            %l = (k^theta_1 * theta_2 * a/W).^(1/(1-theta_2));
            pi = profit(l,a,k);
            inv = investment(k,k_grid);
            %pi = a*k^theta_1*l^theta_2 - W*l;
            phi = ones(1,Nk);
            for ss=1:Nk
                phi(ss) = adjustment_cost(inv(ss),k);
            end
            d = divident(pi,inv,phi);
            %d = pi + (1-delta-b_0)*k-k_grid-b_1*k*(k_grid/k-1).^2;
            [v,ind] = max(d + 1/(1+r)*a_prob(i,:)*value0);
            value(i,j)=v;
            k_policy_index(i,j)=ind;
            k_policy(i,j) = k_grid(k_policy_index(i,j));
            investment_policy4(i,j) = investment(k,k_policy(i,j));
        end
    end
    distance = sum(sum(abs(value-value0)));
    value0=value;
    iteration = iteration + 1;
end
display("iteration =    " + iteration + "   difference =   " + distance)
toc;

%% Plot the whole graph 
optimal_investment_at_mean_shock(4,:)=investment_policy4(m+1,:);

figure(1)
plot(k_grid,optimal_investment_at_mean_shock');
title(['investment policy Nk= ',num2str(Nk)])
ylabel('Optimal Investment at Mean Shock')
xlabel('Current Capital Stock')
legend('low adj. cost','high adj. cost','asymmetric adj. cost','fixed cost');
savefig('q1d_investment_policy') 