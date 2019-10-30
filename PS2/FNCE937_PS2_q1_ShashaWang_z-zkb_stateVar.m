% Shasha Wang PS1 Q3 FNCE 937

% If one wants to do MULTIGRID to speed up the process, just change
% kGridLength to a vector

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
Nb = 30;
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
taxPaymentsFunction = @(k,bond,profit,Rb)ttaoC * (profit - ddelta*k - bond * (Rb-1) * ((1-ttaoC)*profit + (1-ddelta)*k > bond)); % note the non-default indicator
% taxPaymentsFunction = @(k,bond,profit,Rb)ttaoC * (profit - ddelta*k - bond * (Rb-1)); % note the non-default indicator
% taxPaymentsFunction = @(k,bond,profit,Rb)ttaoC * (profit - ddelta*k - bond * (Rb-1)); % without the non-default indicator

dividentFunction = @(profit,investment,bond,bondPrime,RbMinus,taxPayments)(profit - investment  ...
    + bondPrime - RbMinus * bond - taxPayments).*(1 + llambda * ((profit - investment + bondPrime - RbMinus * bond - taxPayments) < 0)); % note the indicator function for issuance cost

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

% Use multigrid method to speed up iteration
kGridLength           = [30]; % number of points in grid for capital
kMin            = 0.000001;
kMax            = 1.5 * kSteadyState;
grid_b = curvspace(kMin,kMax,Nb,5)';

% Required matrices and vectors
% Dimensionality is k,b,a,aMinus
    
kPolicyIndex = zeros(kGridLength(1),Nb,Na,Na);
kPolicy = zeros(kGridLength(1),Nb,Na,Na);
bPolicyIndex = zeros(kGridLength(1),Nb,Na,Na);
bPolicy = zeros(kGridLength(1),Nb,Na,Na);

value   = zeros(kGridLength(1),Nb,Na,Na);
value0 = ones(kGridLength(1),Nb,Na,Na); % initial guess


tic
for i=1:length(kGridLength)
    
    grid_k           = curvspace(kMin,kMax,kGridLength(i),3)';
    
    % Calculate Default Probability
    mCutOffValue = zeros(kGridLength(i),Nb);
    mDefaultProbability = zeros(kGridLength(i),Nb,Na);

    vDenominator = grid_k.^aalphaK .* (aalphaL * grid_k.^aalphaK / W).^(aalphaL/(1-aalphaL)) - W * (aalphaL * grid_k.^aalphaK/W).^(1/(1-aalphaL));

    for ib = 1:Nb
    %     vNumerator = max(0.000000001, (vBond(ib) - (1-ddelta)*k_grid)/(1-ttaoC));
        vNumerator =  (grid_b(ib) - (1-ddelta)*grid_k)/(1-ttaoC);

    %     vCutOffValue = (vNumerator ./v_Denominator).^(1-aalphaL); % I droped
    %     this line because it generates complex numbers

        vCutOffValue = vNumerator ./vDenominator;
        mCutOffValue(:,ib) = vCutOffValue; % cutoff value^(1-aalphaL) is cutoff productivity, but I don't use productivity per se in order to avoid complex numbers
        for ia = 1:Na
            for ik = 1:kGridLength(1)
                 temp= (vCutOffValue(ik) > grid_a.^(1/(1-aalphaL)));
                mDefaultProbability(ik,ib,ia) = sum(temp.* m_a_prob(ia,:)');
            end
        end
    end

    mRf = Rf*ones(kGridLength(1),Nb,Na);
    mCarryOnProbability = 1-mDefaultProbability;
    mRb = min(Rf/mCarryOnProbability,1000000);
    
    tolerance = 0.00001;
    iteration = 0;
    distance = 100;

    kPrime = repmat(grid_k,1,Nb); % Nk*Nb matrix
    bondPrime = repmat(grid_b',kGridLength(i),1); % Nk*Nb matrix

    tic
    while distance > tolerance
%         tic
        for ia=1:Na
            a = grid_a(ia);

            for ik=1:kGridLength(i)
                k = grid_k(ik);        
                labor = laborFunction(a,k);
                profit = profitFunction(a,k,labor);
                investment = investmentFunction(k,kPrime); % Nk*Nb matrix

                for iaMinus = 1:Na
                    RbMinus = mRb(:,:,iaMinus); % Nk*Nb matrix
                    aMinus = grid_a(iaMinus);

                    for ib = 1:Nb
                        bond = grid_b(ib);

                        if (1-ttaoC)*profit + (1-ddelta)*k <= bond
                            value (ik,ib,ia,iaMinus)=0;

                        else

                            taxPayments = taxPaymentsFunction(k,bond,profit,RbMinus); % Nk*Nb matrix
                            divident = dividentFunction(profit,investment,bond,bondPrime,RbMinus,taxPayments); % Nk*Nb matrix

                            valueTomorrow = zeros(kGridLength(i),Nb,Na);
                            for iaPrime = 1:Na % iterate over all possible states for tomorrow
                                aPrime = grid_a(iaPrime);
                                laborPrime = laborFunction(aPrime,kPrime);
                                profitPrime = profitFunction(aPrime,kPrime,laborPrime);
                                valueTomorrow(:,:,iaPrime) = value0(:,:,iaPrime,ia) * m_a_prob(ia,iaPrime) .* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime > bondPrime);% 需要考虑default之后value为0
                            end
                            valueTomorrow = sum(valueTomorrow,3); % sum by the third dimension to get a Nk*Nb matrix

                            x = divident + M * valueTomorrow;

                            [rows,cols]=find(x==max(max(x)));

    %                 if (1-ttaoC)*profit + (1-ddelta)*k <= bond
    %                     value(ik,ib,ia) = 0;
    %                 else
    %                     value(ik,ib,ia) = x(max(rows),max(cols));
    %                 end

                            kPolicyIndex(ik,ib,ia,iaMinus) = max(rows);
                            bPolicyIndex(ik,ib,ia,iaMinus) = max(cols);

                            kPolicy(ik,ib,ia,iaMinus) = grid_k(max(rows));
                            bPolicy(ik,ib,ia,iaMinus) = grid_b(max(cols));
                            value(ik,ib,ia,iaMinus) = max(max(x));
                        end
                    end
                end
            end
        end

%         toc
        %hold on
        %mesh(aa, kk, value')
        distance = sum(sum(sum(sum(abs(value(:,:,:,:)-value0(:,:,:,:))))));
        value0 = value;
        iteration = iteration + 1;

        if mod(iteration,5) == 0
            display("iteration =    " + iteration + "   difference =   " + distance )
        end
    end

    display("iteration =    " + iteration + "   difference =   " + distance + ". Converged")
    if i ~= length(kGridLength)
        value0 = interp1(grid_k,value,linspace(kMin, kMax, kGridLength(i+1)));
        value  = value0;
%         kPolicy = interp1(grid_k,kPolicy,linspace(kMin, kMax, kGridLength(i+1)));
        kPolicy         = zeros(kGridLength(i+1),Nb,Na,Na);
        kPolicyIndex    = zeros(kGridLength(i+1),Nb,Na,Na);
        bPolicy         = zeros(kGridLength(i+1),Nb,Na,Na);
        bPolicyIndex    = zeros(kGridLength(i+1),Nb,Na,Na);
        
    end
    
end

% I left the laptop in the office to let it run.
%    "iteration =    1636   difference =   9.9406e-06. Converged"

toc
save('valuePrevious','value')

figure(3);
[bb,kk]=meshgrid(grid_b, grid_k);
mesh(bb, kk, value(:,:,1,round((Na+1)/2)));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, value(:,:,ia,round((Na+1)/2)));    
end

title('Value Under Different Shocks given mean $z^{-}$','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('Value','interpreter','latex')
savefig('q1c_value_3D')

figure(4)
mesh(bb, kk, kPolicy(:,:,1,round((Na+1)/2)));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, kPolicy(:,:,ia,round((Na+1)/2)));    
end

title('Policy $k^\prime$ Under Different Shocks given mean $z^{-}$','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('$k^\prime$','interpreter','latex')
savefig('q1c_kPolicy_3D')


figure(5)
mesh(bb, kk, bPolicy(:,:,1,round((Na+1)/2)));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, bPolicy(:,:,ia,round((Na+1)/2)));    
end

title('Policy $b^\prime$ Under Different Shocks given mean $z^{-}$','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('$bond^\prime$','interpreter','latex')
savefig('q1c_bPolicy_3D')

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





    
    



