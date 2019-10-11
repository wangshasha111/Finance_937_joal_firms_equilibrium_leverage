% Shasha Wang PS1 Q3 FNCE 937
clear all;
close all;

cd 'E:\Dropbox\fall 19-20\Finance 937\PS1\question_3'

%% Q3 Corporate Investment in Equilibrium - Gomes(2001)

%% Parametization
M = 0.95;
r = 1/M-1; % interest rate notation convenience
delta = 0.1;
alpha_k = 0.3;
alpha_l = 0.65;
W = 2; % wage
rho = 0.95;
sigma = 0.02;
m = 4;
%Na = 2*m+1;
Na = 15;


% exit rate calibrated to 2.5%

%% k_grid
Nk = 200;
k_min = 0.0001;
k_SteadyState = 1;
a_mean = ((r+delta)/alpha_k)^(1-alpha_l)*(W/alpha_l)^alpha_l; % set z such that k_SteadyState=1;
log_a_mean = log(a_mean);
log_a_min=log(a_mean)-m*sigma/sqrt(1-rho^2);
log_a_max=log(a_mean)+m*sigma/sqrt(1-rho^2);
log_a_grid=linspace(log_a_min,log_a_max,Na);
a_min=exp(log_a_min);
a_max=exp(log_a_max);
k_max = (a_max/delta)^(1/(1-alpha_k));
k_max = min(2*k_SteadyState, k_max); % Tighten the grid
%k_grid = linspace(k_min, k_max, Nk);
k_grid = curvspace(k_min,k_max,Nk,2); % I use curved grid to enhance accuracy

%% a's transition matrix

[A,a_prob] = tauchen(Na,log_a_mean,rho,sigma,m);
assert(sum(abs(A-log_a_grid))<0.01);

a_grid = exp(log_a_grid);

clear A;

%% create functions for convenience
labor_function = @(a,k) (k.^alpha_k * alpha_l * a/W).^(1/(1-alpha_l));
profit_function = @(a,k,labor,fixed_cost)a * k.^alpha_k.* labor.^alpha_l - W*labor - fixed_cost;
investment_function = @(k,k_prime)k_prime - (1-delta)*k; %k_prime usually is k_grid
adjustment_cost_function = @(investment,k)0.5*k*(investment/k-delta).^2;
divident_function = @(profit,investment,adjustment_cost)profit-investment-adjustment_cost;

%% Value Function Iteration and Policy Function



%%
Nfixed_cost = 10;
% fixed_cost_min = 0.050301111111111;
% fixed_cost_max = 0.050346666666667;
%fixed_cost_min = 0.02393;
%fixed_cost_max = 0.02397;
fixed_cost_min = 0.02337;
fixed_cost_max = 0.02339;
fixed_cost_grid = linspace(fixed_cost_min,fixed_cost_max,Nfixed_cost);
entry_rate_conditional_grid = zeros(1,Nfixed_cost);
exit_rate_conditional_grid = zeros(1,Nfixed_cost);
entry_rate_unconditional_grid = zeros(1,Nfixed_cost);
exit_rate_unconditional_grid = zeros(1,Nfixed_cost);
exit_rate_desired = 0.025;

tic;
for n=1:Nfixed_cost
    k_policy_index = zeros(Na,Nk);
    k_policy = zeros(Na,Nk);
    investment_policy = zeros(Na,Nk);
    value   = zeros(Na,Nk);
    value0 = ones(Na,Nk); % initial guess
    tolerance = 0.00001;
    iteration = 0;
    distance = 100;

    [aa,kk]=meshgrid(a_grid, k_grid);

    %figure(1)
    %mesh(aa, kk, value0')

    while distance > tolerance
        for i=1:Na
            for j=1:Nk
                a = a_grid(i);
                k = k_grid(j);
                labor = labor_function(a,k);
                profit = profit_function(a,k,labor,fixed_cost_grid(n));
                investment = investment_function(k,k_grid);
                adjustment_cost = adjustment_cost_function(investment,k);
                divident = divident_function(profit,investment,adjustment_cost);
                [v,ind] = max(divident + M*max(a_prob(i,:)*value0,0));
                value(i,j)=v;
                k_policy_index(i,j)=ind;
                k_policy(i,j) = k_grid(k_policy_index(i,j));
            end
        end
        %hold on
        %mesh(aa, kk, value')
        distance = sum(sum(abs(value-value0)));
        value0=value;
        iteration = iteration + 1;
    end
    display("iteration =    " + iteration + "   difference =   " + distance + "   fixed cost =   " + fixed_cost_grid(n))

    %% transition prob on z*k-by-z*k space
    transition_matrix = zeros(Na*Nk,Na*Nk);
    for i=1:Na
        k_prob = zeros(Nk,Nk);
        for j=1:Nk
            j1 = k_policy_index(i,j);
            k_prob(j,j1)=1;
        end
        transition_matrix((i-1)*Nk+1:i*Nk,:) = kron(a_prob(i,:),k_prob);
    end

    %% calculate invariant transition
    distribution0=( 1/(Nk*Na) )*ones(1,Nk*Na); % initial guess
    distance=100; tolerance=0.0001;
    iteration = 0;
    while distance>tolerance
        distribution = distribution0*transition_matrix;
        distance=sum(abs(distribution-distribution0));
        distribution0 = distribution;
        iteration = iteration+1;
    end

    % transform distribution vector into a matrix
    distribution_Na_By_Nk=reshape(distribution,[Na,Nk]);

    % compute entry rate by looking at potential entrants
    entry_rate_conditional_grid(1,n) = sum((k_policy_index(:,1)>1).* (distribution_Na_By_Nk(:,1)))/sum(distribution_Na_By_Nk(:,1));
    exit_rate_conditional_grid(1,n) = sum(sum((k_policy_index(:,2:Nk)==1).* (distribution_Na_By_Nk(:,2:Nk))))/sum(sum(distribution_Na_By_Nk(:,2:Nk)));
    entry_rate_unconditional_grid(1,n) = sum((k_policy_index(:,1)>1).* (distribution_Na_By_Nk(:,1)));
    exit_rate_unconditional_grid(1,n) = sum(sum((k_policy_index(:,2:Nk)==1).* (distribution_Na_By_Nk(:,2:Nk))));

end
toc;


close all;
figure(1);
subplot(2,1,1)
plot(fixed_cost_grid,exit_rate_conditional_grid);
yline(exit_rate_desired);
%ylim([0,max(exit_rate_grid)])
xlim([fixed_cost_min,fixed_cost_max])
xlabel('fixed cost');
ylabel('conditional exit rate');
title('conditional exit rate under different fixed cost');
%savefig('q3a_fixed_cost_exit_rate_conditional')

subplot(2,1,2)
plot(fixed_cost_grid,exit_rate_unconditional_grid);
yline(exit_rate_desired);
%ylim([0,max(exit_rate_grid)])
xlim([fixed_cost_min,fixed_cost_max])
xlabel('fixed cost');
ylabel('unconditional exit rate');
title('unconditional exit rate under different fixed cost');
savefig('q3a_fixed_cost_exit_rate')

[v,ind]=min(abs(exit_rate_unconditional_grid-exit_rate_desired));

%% Plot the stationary distribution when the measure of firms is 1
figure(2)
[aa,kk]=meshgrid(a_grid, k_grid);
mesh(aa, kk, distribution_Na_By_Nk')
title('firm stationary distribution')
xlabel('Productivity')
ylabel('Current Capital Stock')
zlabel('Probability Mass')
zlim([0,max(max(distribution_Na_By_Nk))]);
savefig('q3a_stationary_distribution') 

%% Compute Aggregate Labor demand
% labor_function = @(a,k) (k.^alpha_k * alpha_l * a/W).^(1/(1-alpha_l));
labor_demand = zeros(Na,Nk);
for ai=1:Na
    for ki=1:Nk
        labor_demand(ai,ki)=labor_function(a_grid(ai),k_grid(ki));
    end
end

aggregate_labor_demand = sum(sum(distribution_Na_By_Nk .* labor_demand));

% in equilibrium supply=demand
B = aggregate_labor_demand/W^0.1;

% entry rate
fixed_cost=fixed_cost_grid(ind);
entry_rate=entry_rate_unconditional_grid(ind);
exit_rate=exit_rate_unconditional_grid(ind);

display("Fixed cost =    " + fixed_cost );
display("Exit rate =    " + exit_rate );
display("Entry rate =    " + entry_rate );
display("B =    " + B );

%% cost of entry
% compute the invariant distribution of a_grid
distribution_a0=( 1/Na )*ones(1,Na); % initial guess
distribution_a=distribution_a0;
distance=100; tolerance=0.00001;
iteration = 0;
while distance>tolerance
    distribution = distribution_a*a_prob;
    distance=sum(abs(distribution-distribution_a));
    distribution_a = distribution;
    iteration = iteration+1;
end

% conditional on (a_i,k),compute the expected value of state (a_j,k(a_i,k))
value_next_period = zeros(1,Na);
value_entry = zeros(1,Na);
for ai = 1:Na % this period shock
    k_prime_index = k_policy_index(ai,1);
    for aii = 1:Na % next period shock
        value_next_period(1,aii) = value(aii,k_prime_index);
    end
    value_entry(1:ai) = sum(a_prob(ai,:).* value_next_period); % conditional on a_i, the expected value
end

% then compute the expected value
expected_value_entry = sum(value_entry .* distribution_a);

display('Setting next period capital stock according to this period shock, ')
display("Entry cost =    " + expected_value_entry  );
assert(expected_value_entry>0);
display('As expected, the entry cost is positive.');

%% Using LG's method of computing entry cost
% conditional on (a_i,k),compute the expected value of state (a_j,k(a_i,k))
value_next_period = zeros(1,Na);
value_entry = zeros(1,Na);
for ai = 1:Na % this period shock
    for aii = 1:Na % next period shock
        k_prime_index = k_policy_index(aii,1);
        value_next_period(1,aii) = value(aii,k_prime_index);
    end
    value_entry(1:ai) = sum(a_prob(ai,:).* value_next_period); % conditional on a_i, the expected value
end

% then compute the expected value
expected_value_entry_LG = sum(value_entry .* distribution_a);

display('Using another method suggested by LG, which set next period capital stock according to the next period shock,')
display("Entry cost =    " + expected_value_entry_LG  );
assert(expected_value_entry_LG>0);
display('As expected, the entry cost is positive.');

assert(abs(expected_value_entry_LG-expected_value_entry)<0.00001);
display('Please note that the results are the same. The two methods are not necessairily equivalent, but the difference in results is too subtle to be significant.');


%% q3b suppose B=1; What is the mass of firms
B2 = 1;
aggregate_labor_supply = B2 * W^0.1;
mass_firm = aggregate_labor_supply/aggregate_labor_demand;
display("Mass of firms =    " + mass_firm  );
