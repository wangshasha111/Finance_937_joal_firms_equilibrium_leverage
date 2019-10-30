% Shasha Wang PS1 Q3 FNCE 937
% Na = 15, better set to be this, otherwise the legend of figure 2 would
% adjust

clear;
close all;

cd 'E:\Dropbox\fall 19-20\Finance 937\PS2\question_1'

%% Part I: Liquidity/Covenant Default

% Parametization
M = 0.99;
Rf = 1/M;
r = 1/M-1; % interest rate for notation convenience
ddelta = 0.1;
aalphaK = 0.3;
aalphaL = 0.6;
W = 2; % wage
llambda = 0.025; % proportional cost of issuing equity
ttaoC=0.15; % corporate tax rate

%% k_grid & b_grid & a_grid & m_a_prob （a's transition probability matrix）
Nk = 20;
kMin = 0.0001;

kSteadyState = 1;
tempK = kSteadyState^((aalphaK+aalphaL-1)/(1-aalphaL));
aMean = ((r+ddelta)/aalphaK/tempK)^(1-aalphaL)*(W/aalphaL)^aalphaL; % set a such that k_SteadyState equals how much you set it to be;

% a_grid and a's transition matrix
m = 3; % parameter for tauchen method
Na = 3;
rrho = 0.7;
ssigma = 0.05;
[grid_a_log,m_a_prob] = tauchen(Na,log(aMean),rrho,ssigma,m);
grid_a = exp(grid_a_log)';
grid_a_minus = grid_a;
aMax = max(grid_a);

% k_grid
kMax = (aMax/ddelta)^(1/(1-aalphaK));
kMax = min(2*kSteadyState, kMax); % Tighten the grid
grid_k = curvspace(kMin,kMax,Nk,2)'; % I use curved grid to enhance accuracy

% b_grid grid for bond
% b_grid should be finer to see the difference in default probability under different productivity shocks at steady state
Nb = 15;
grid_b = curvspace(kMin,kMax,Nb,5)'; % To cover up as wide leverage level as possible

% invariant distribution of a_grid - vDistribution_a0
vDistribution_a0=( 1/Na )*ones(Na,1); % initial guess
vDistribution_a=vDistribution_a0;
distance=100; tolerance=0.00001;
iteration = 0;
while distance>tolerance
    distribution = vDistribution_a'*m_a_prob;
    distance=sum(abs(distribution-vDistribution_a'));
    vDistribution_a = distribution';
    iteration = iteration+1;
end

clear vDistribution_a0 distribution;

%% Functions
% create functions for convenience
laborFunction = @(a,k) (k.^aalphaK * aalphaL * a/W).^(1/(1-aalphaL));
profitFunction = @(a,k,labor)a * k.^aalphaK.* labor.^aalphaL - W*labor;
% profit = a^(1/(1-aalphaL)) * (k.^aalphaK * (aalphaL * k.^aalphaK / W).^(aalphaL/(1-aalphaL)) - W * (aalphaL * k.^aalphaK/W).^(1/(1-aalphaL)));
investmentFunction = @(k,kPrime)kPrime - (1-ddelta)*k; %k_prime usually is k_grid
taxPaymentsFunction = @(k,bond,profit,Rb)ttaoC * (profit - ddelta*k - bond * (Rb-1) * ((1-ttaoC)*profit + (1-ddelta)*k >= bond)); % note the non-default indicator
% taxPaymentsFunction = @(k,bond,profit,Rb)ttaoC * (profit - ddelta*k - bond * (Rb-1)); % without the non-default indicator

dividentFunction = @(profit,investment,bond,bondPrime,Rb,taxPayments)(profit - investment  ...
    + bondPrime - Rb * bond - taxPayments).*(1 + llambda * ((profit - investment + bondPrime - Rb * bond - taxPayments) < 0)); % note the indicator function for issuance cost

%% a) Probability of Default Next Period
% Compute and plot the probability of default next period, ...
% conditional on the value of the shocks today p(z,b',k')

% rewrite profit = a^(1/(1-aalphaL)) * (k.^aalphaK .* (aalphaL * k.^aalphaK / W).^(aalphaL/(1-aalphaL)) - W * (aalphaL * k.^aalphaK/W).^(1/(1-aalphaL)));
% once we plug in the FOC for labor

% vLabor = laborFunction(a,k_grid);
% vProfit = profitFunction(a,k_grid,vLabor);

mCutOffValue = zeros(Nk,Nb);
mDefaultProbability = zeros(Nk,Nb,Na);

vDenominator = grid_k.^aalphaK .* (aalphaL * grid_k.^aalphaK / W).^(aalphaL/(1-aalphaL)) - W * (aalphaL * grid_k.^aalphaK/W).^(1/(1-aalphaL));

for ib = 1:Nb
%     vNumerator = max(0.000000001, (vBond(ib) - (1-ddelta)*k_grid)/(1-ttaoC));
    vNumerator =  (grid_b(ib) - (1-ddelta)*grid_k)/(1-ttaoC);

%     vCutOffValue = (vNumerator ./v_Denominator).^(1-aalphaL); % I droped
%     this line because it generates complex numbers

    vCutOffValue = vNumerator ./vDenominator;
    mCutOffValue(:,ib) = vCutOffValue; % cutoff value^(1-aalphaL) is cutoff productivity, but I don't use productivity per se in order to avoid complex numbers
    for ia = 1:Na
        for ik = 1:Nk
             temp= (vCutOffValue(ik) > grid_a.^(1/(1-aalphaL)));
            mDefaultProbability(ik,ib,ia) = sum(temp.* m_a_prob(ia,:)');
        end
    end
end
            
% plot the 3D matrix layer by layer

figure(1);
[bb,kk]=meshgrid(grid_b, grid_k);
mesh(bb, kk, mDefaultProbability(:,:,1));

for ia = 2:Na
    hold on;
    mesh(bb, kk, mDefaultProbability(:,:,ia));    
end

title('Next Period Default Probability')
ylabel('Next Period Capital Stock $k^\prime$','interpreter','latex')
xlabel('Next Period Debt $b^\prime$','interpreter','latex')
zlabel('Next Period Default Probability','interpreter','latex')
zlim([0,1])
savefig('q1a_default_prob_3D')

% plot a 2D default risk graph in steady state k=1
% we cannot tell from the graph differences among different productivity
% situations
figure(2)
[v,ind] = max(-abs(grid_k - kSteadyState));
plot(grid_b,mDefaultProbability(ind,:,1));
% for ia=[(Na+1)/2,Na]
for ia=2:Na
    hold on
    plot(grid_b,mDefaultProbability(ind,:,ia));
end
ylim([0,1]);
% xlim([0.5,2]);
% xlim([0.75,1.2]);
legend('low productivity','median productivity','high productivity','Location','southeast')

title('Default Probability at Steady State $p(z,k^\prime,b^\prime)$','interpreter','latex')
xlabel('Next Period Debt $b^\prime$','interpreter','latex')
ylabel('Next Period Default Probability','interpreter','latex')
savefig('q1a_default_prob_2D')

%% b) Required Rate of Return of Bonds
% Use the conditional default probability to compute the required rate of return ...
% by bondholders, Rb(b'; k'; z) that ensures they make 0 profits. ...
% Assume for simplicity bondholders get paid 0 upon default.

mRf = Rf*ones(Nk,Nb,Na);
mCarryOnProbability = 1-mDefaultProbability;
mRb = min(Rf/mCarryOnProbability,1000000);
% mRb = mRf./mCarryOnProbability; % cross this out because we have to deal
% with denominator being zero

% for ia = 1:Na
%     for ib = 1:Nb
%         for ik = 1:Nk
%             if mCarryOnProbability(ik,ib,ia)==0
%                 mRb(ik,ib,ia) =10000000000; % deal with complex numbers
%             else
%                 mRb(ik,ib,ia) = Rf/mCarryOnProbability(ik,ib,ia);
%             end
%         end
%     end
% end

%% c) Value Function Iteration
% Solve the Bellman equation for the equity holders taking as given the
% function for the required rate of return by bondholders, Rb

% mMixKandB=zeros(Nk,Nb);
% for ib = 1:Nb
%     for ik = 1:Nk
%         mMixKandB(ik,ib)=[grid_k(ik),grid_b(ib)];
%     end
% end

kPolicyIndex = zeros(Nk,Nb,Na,Na);
kPolicy = zeros(Nk,Nb,Na,Na);
bPolicyIndex = zeros(Nk,Nb,Na,Na);
bPolicy = zeros(Nk,Nb,Na,Na);

value   = zeros(Nk,Nb,Na,Na);
% value0 = repmat([1:Nb],Nk,1,Na);
value0 = ones(Nk,Nb,Na,Na); % initial guess

    %[aa,kk]=meshgrid(a_grid, k_grid);

    %figure(1)
    %mesh(aa, kk, value0')
    
tolerance = 0.00001;
iteration = 0;
distance = 100;

kPrime = repmat(grid_k,1,Nb); % Nk*Nb matrix
bondPrime = repmat(grid_b',Nk,1); % Nk*Nb matrix

tic
while distance > tolerance
    tic
    for ia=1:Na
        Rb = mRb(:,:,ia); % Nk*Nb matrix
        a = grid_a(ia);
        for ib = 1:Nb
            for ik=1:Nk
                k = grid_k(ik);
                bond = grid_b(ib);
                
                labor = laborFunction(a,k);
                profit = profitFunction(a,k,labor);

                investment = investmentFunction(k,kPrime); % Nk*Nb matrix
%                 investment = repmat(investment,1,Nb); % Nk*Nb matrix
%                 bondPrime = repmat(grid_b',Nk,1);
%                 adjustment_cost = adjustment_cost_function(investment,k);
%                 tax = taxPaymentsFunction(a,k,bond,profit,mRb(ik,ib,ia));
                taxPayments = taxPaymentsFunction(k,bond,profit,Rb); % Nk*Nb matrix
                divident = dividentFunction(profit,investment,bond,bondPrime,Rb,taxPayments); % Nk*Nb matrix
               
                valueTomorrow = zeros(Nk,Nb,Na);
                for iaa = 1:Na
                    valueTomorrow(:,:,iaa) = value0(:,:,iaa) * m_a_prob(ia,iaa);% 这里有问题，需要考虑default之后value为0
                end
                valueTomorrow = sum(valueTomorrow,3); % sum by the third dimension
                
                x = divident + M * valueTomorrow;
                
                [rows,cols]=find(x==max(max(x)));
                
%                 if (1-ttaoC)*profit + (1-ddelta)*k <= bond
%                     value(ik,ib,ia) = 0;
%                 else
%                     value(ik,ib,ia) = x(max(rows),max(cols));
%                 end

                kPolicyIndex(ik,ib,ia) = max(rows);
                bPolicyIndex(ik,ib,ia) = max(cols);
                
                kPolicy(ik,ib,ia) = grid_k(max(rows));
                bPolicy(ik,ib,ia) = grid_b(max(cols));
                
            end
        end
    end
    toc
    %hold on
    %mesh(aa, kk, value')
    distance = sum(sum(sum(abs(value(:,:,:)-value0(:,:,:)))));
    value0 = value;
    iteration = iteration + 1;

    if mod(iteration,5) == 0
        display("iteration =    " + iteration + "   difference =   " + distance )
    end
end

display("iteration =    " + iteration + "   difference =   " + distance + ". Converged")

toc

figure(3);
[bb,kk]=meshgrid(grid_b, grid_k);
mesh(bb, kk, value(:,:,1));

for ia = 2:Na
    hold on;
    mesh(bb, kk, value(:,:,ia));    
end

title('Value Under Different Shocks')
ylabel('Capital Stock $k^\prime$','interpreter','latex')
xlabel('Debt $b^\prime$','interpreter','latex')
zlabel('Value','interpreter','latex')
savefig('q1c_value_3D')

%% d) Stationary Distribution of Firms
% Consider now a world with many such firms and no entry or exit. Specifically,
% suppose that upon hitting the default threshold debt claims are settled so b = 0.
% The restructured firm continues to operate but with capital, k = 0 and the previous
% productivity shock, z.

% Compute the stationary distribution of firms. Use this distribution to construct
% a table reporting the cross-sectional average values of:
% (1) probability of default, p(·);
% (2) required return on risky bonds, Rb(·);
% (3) leverage ratio, b=k;
% (4) investment to capital ratio, i=k.
% (5) fraction of firms issuing equity;





    
    



