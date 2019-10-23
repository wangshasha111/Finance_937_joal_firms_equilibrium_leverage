% Shasha Wang PS1 Q3 FNCE 937
clear;
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

New_Entrant=20;

% exit rate calibrated to 2.5%

%% k_grid & a_grid
Nk = 100;
k_min = 0.0001;
% k_max=1.5;
% k_min = 0;
% k_SteadyState = 1;
% a_mean = ((r+delta)/alpha_k)^(1-alpha_l)*(W/alpha_l)^alpha_l; % set a such that k_SteadyState=1;

k_SteadyState = 1;
temp_k = k_SteadyState^((alpha_k+alpha_l-1)/(1-alpha_l));
a_mean = ((r+delta)/alpha_k/temp_k)^(1-alpha_l)*(W/alpha_l)^alpha_l; % set a such that k_SteadyState=2;

log_a_mean = log(a_mean);
log_a_min=log(a_mean)-m*sigma/sqrt(1-rho^2);
log_a_max=log(a_mean)+m*sigma/sqrt(1-rho^2);
log_a_grid=linspace(log_a_min,log_a_max,Na);
a_min=exp(log_a_min);
a_max=exp(log_a_max);
k_max = (a_max/delta)^(1/(1-alpha_k));
k_max = min(11*k_SteadyState, k_max); % Tighten the grid
%k_grid = linspace(k_min, k_max, Nk);
k_grid = curvspace(k_min,k_max,Nk,2); % I use curved grid to enhance accuracy

%% a's transition matrix

[A,a_prob] = tauchen(Na,log_a_mean,rho,sigma,m);
assert(sum(abs(A-log_a_grid))<0.01);

a_grid = exp(log_a_grid);

clear A;

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

clear distribution_a0 distribution;

%% entry distribution gamma implied by the invariant distribution of a
temp = [distribution_a;zeros(Nk-1,Na)];
entry_distribution = reshape(temp,1,Na*Nk);% reshape temp into a Na*Nk row vector

%% create functions for convenience
labor_function = @(a,k) (k.^alpha_k * alpha_l * a/W).^(1/(1-alpha_l));
profit_function = @(a,k,labor,fixed_cost)a * k.^alpha_k.* labor.^alpha_l - W*labor - fixed_cost;
investment_function = @(k,k_prime)k_prime - (1-delta)*k; %k_prime usually is k_grid
adjustment_cost_function = @(investment,k)0.5*k*(investment/k-delta).^2;
divident_function = @(profit,investment,adjustment_cost)profit-investment-adjustment_cost;

%% Value Function Iteration and Policy Function
% loop over different fixed cost

Nfixed_cost = 10;
% fixed_cost_min = 0.050301111111111;
% fixed_cost_max = 0.050346666666667;
% fixed_cost_min = 0.02393;
% fixed_cost_max = 0.02397;

% fixed_cost_min = 0.02337;
% fixed_cost_max = 0.02339;

% fixed_cost_min = 0.18;
% fixed_cost_max = 0.2222;

% fixed_cost_min = 0.17;
% fixed_cost_max = 0.1847;

% fixed_cost_min = 0.1798;
% fixed_cost_max = 0.1814;

% fixed_cost_min = 0.8;
% fixed_cost_max = 1.0;

% fixed_cost_min = 0.05;
% fixed_cost_max = 0.15;

% fixed_cost_min = 0.09691;
% fixed_cost_max = 0.1;

fixed_cost_min = 0.025;
fixed_cost_max = 0.03;

fixed_cost_grid = linspace(fixed_cost_min,fixed_cost_max,Nfixed_cost);
%entry_rate_conditional_grid = zeros(1,Nfixed_cost);
%exit_rate_conditional_grid = zeros(1,Nfixed_cost);
%entry_rate_unconditional_grid = zeros(1,Nfixed_cost);
%exit_rate_unconditional_grid = zeros(1,Nfixed_cost);

entry_rate_grid = zeros(1,Nfixed_cost);
exit_rate_grid = zeros(1,Nfixed_cost);
exit_rate_desired = 0.025;
distribution_Na_By_Nk=zeros(Na,Nk,Nfixed_cost); % a 3D matrix
k_policy_index = zeros(Na,Nk,Nfixed_cost);
k_policy = zeros(Na,Nk,Nfixed_cost);
exit_policy = zeros(Na,Nk,Nfixed_cost);% if firm chooses to exit, entry is set to 1
stay_policy = zeros(Na,Nk,Nfixed_cost);% if firm chooses to stay, entry is set to 1
investment_policy = zeros(Na,Nk,Nfixed_cost);
value   = zeros(Na,Nk,Nfixed_cost);
value0 = ones(Na,Nk,Nfixed_cost); % initial guess
transition_matrix = zeros(Na*Nk,Na*Nk,Nfixed_cost);
%distribution0=( 1/(Nk*Na) )*ones(1,Nk*Na,Nfixed_cost); % initial guess of invariant distribution
distribution = zeros(1,Nk*Na,Nfixed_cost); % initial guess of invariant distribution


tic;
for nn=1:Nfixed_cost

    tolerance = 0.00001;
    iteration = 0;
    distance = 100;

    %[aa,kk]=meshgrid(a_grid, k_grid);

    %figure(1)
    %mesh(aa, kk, value0')

    while distance > tolerance
        for j=1:Nk
            for i=1:Na
                a = a_grid(i);
                k = k_grid(j);
                labor = labor_function(a,k);
                profit = profit_function(a,k,labor,fixed_cost_grid(nn));
                investment = investment_function(k,k_grid);
                adjustment_cost = adjustment_cost_function(investment,k);
                divident = divident_function(profit,investment,adjustment_cost);
                [v,ind] = max(divident + M*max(a_prob(i,:)*value0(:,:,nn), 0));
                value(i,j,nn)=v;              
                k_policy_index(i,j,nn)=ind;                
                k_policy(i,j,nn) = k_grid(k_policy_index(i,j,nn));
                exit_policy(i,j,nn) = (a_prob(i,:)*value0(:,ind,nn) <= 0);
                stay_policy(i,j,nn) = (a_prob(i,:)*value0(:,ind,nn) > 0);
            end
        end
        %hold on
        %mesh(aa, kk, value')
        distance = sum(sum(abs(value(:,:,nn)-value0(:,:,nn))));
        value0(:,:,nn)=value(:,:,nn);
        iteration = iteration + 1;
    end
    display("loop  " + nn);
    display("iteration =    " + iteration + "   difference =   " + distance + "   fixed cost =   " + fixed_cost_grid(nn))
    
    %% transition prob on z*k-by-z*k space
    for i=1:Na
        k_prob = zeros(Nk,Nk);
        for j=1:Nk
            j1 = k_policy_index(i,j,nn);
            if stay_policy(i,j,nn) == 1
                k_prob(j,j1)=1;
            else
                k_prob(j,j1)=0;
            end
        end
        transition_matrix((i-1)*Nk+1:i*Nk,:,nn) = kron(a_prob(i,:),k_prob);
    end
    
    
    %% calculate invariant transition
%      distribution(:,:,nn) = inv(speye(Na*Nk) - transition_matrix(:,:,nn).*stay_policy(:,:,nn))* entry_distribution' * New_Entrant;  
%     distribution(:,:,nn) = inv(speye(Na*Nk) - transition_matrix(:,:,nn))* entry_distribution' * New_Entrant;  
     distribution(:,:,nn) = entry_distribution'\inv(speye(Na*Nk) - transition_matrix(:,:,nn))* New_Entrant;  
     
     % scale the distribution to measure 1
     distribution(:,:,nn) = distribution(:,:,nn)/sum(distribution(:,:,nn));
     
     distribution_Na_By_Nk(:,:,nn)=reshape(distribution(:,:,nn),[Nk,Na])';
    
    exit_rate_grid(1,nn) = sum(sum(exit_policy(:,:,nn) .* (distribution_Na_By_Nk(:,:,nn))));
    display("exit rate =  " + exit_rate_grid(1,nn));

end
toc;

%% plot
close all;

figure(1);
plot(fixed_cost_grid,exit_rate_grid);
yline(exit_rate_desired);
%ylim([0,max(exit_rate_grid)])
xlim([fixed_cost_min,fixed_cost_max])
xlabel('fixed cost');
ylabel('exit rate');
title('exit rate under different fixed cost');
savefig('q3a_fixed_cost_exit_rate_v1')

% [v,ind]=min(abs(exit_rate_unconditional_grid-exit_rate_desired));
[v,ind]=min(abs(exit_rate_grid-exit_rate_desired));

%% q3a Plot the stationary distribution when the measure of firms is 1
figure(2)
[aa,kk]=meshgrid(a_grid, k_grid);
mesh(aa, kk, distribution_Na_By_Nk(:,:,ind)')
title('firm stationary distribution')
xlabel('Productivity')
ylabel('Current Capital Stock')
zlabel('Probability Mass')
zlim([min(min(min(distribution_Na_By_Nk(:,:,ind)))),max(max(max(distribution_Na_By_Nk(:,:,ind))))])

%zlim([0,max(max(distribution_Na_By_Nk(:,:,ind)))]);
savefig('q3a_stationary_distribution_v1') 

%% Compute Aggregate Labor demand
% labor_function = @(a,k) (k.^alpha_k * alpha_l * a/W).^(1/(1-alpha_l));
labor_demand = zeros(Na,Nk);
for ai=1:Na
    for ki=1:Nk
        labor_demand(ai,ki)=labor_function(a_grid(ai),k_grid(ki));
    end
end

aggregate_labor_demand = sum(sum(distribution_Na_By_Nk(:,:,ind) .* labor_demand));

% In equilibrium, supply==demand
B = aggregate_labor_demand/W^0.1;

% entry rate
fixed_cost=fixed_cost_grid(ind);
%entry_rate=entry_rate_unconditional_grid(ind);
%exit_rate=exit_rate_unconditional_grid(ind);
exit_rate=exit_rate_grid(ind);

display("Fixed cost =    " + fixed_cost );
display("Exit rate =    " + exit_rate );
% display("Entry rate =    " + entry_rate );
display("B =    " + B );

%% cost of entry

% compute the expected value
expected_value_entry = sum(value(:,1,ind)' .* distribution_a);
entry_cost = expected_value_entry;

display('Setting next period capital stock according to this period shock, ')
display("Entry cost =    " + entry_cost  );
assert(entry_cost>0);
display('As expected, the entry cost is positive.');

% compute entry rate
entry_decision = (value(:,1,ind)>entry_cost);
entry_rate = sum(entry_decision'.*distribution_a);

display("Entry rate =    " + entry_rate );

%% q3b suppose B=1; What is the mass of firms
B2 = 1;
aggregate_labor_supply = B2 * W^0.1;
mass_firm = aggregate_labor_supply/aggregate_labor_demand;
display("Mass of firms =    " + mass_firm  );
