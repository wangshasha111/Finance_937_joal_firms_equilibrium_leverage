% Shasha Wang PS2 Q3 FNCE 937

% If one wants to do MULTIGRID to speed up the process, just change
% kGridLength to a vector

clear;
close all;

% cd 'E:\Dropbox\fall 19-20\Finance 937\PS2\question_3'

%% Part III: Optimal Default
% First we need to decide on the default behavior - once default, will the
% firm continue to run the firm as the ownership with 0 capital and 0 bond 
% as in question 1(d), or will its ownership be changed and thus value put 
% to 0? 

% It is implied in the Bellman equation that, once default, firms exit the
% market, since next period we are choosing between 0 and e(z',k',b'), not
% between e(z',0,0) and e(z',k',b').

% Parametization
M = 0.99;
Rf = 1/M;
RbConstant = 1.01 * Rf;
r = 1/M-1; % risk free interest rate for notation convenience
ddelta = 0.1;
aalphaK = 0.3;
aalphaL = 0.6;
W = 2; % wage
llambda = 0.025; % proportional cost of issuing equity
ttaoC=0.15; % corporate tax rate

%% a_grid & m_a_prob （a's transition probability matrix）
% Nk = 20;
% kMin = 0.00001;

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
% grid_a_minus = grid_a;
aMax = max(grid_a);

% k_grid
kMax = (aMax/ddelta)^(1/(1-aalphaK));
kMax = min(2*kSteadyState, kMax); % Tighten the grid
% grid_k = curvspace(kMin,kMax,Nk,2)'; % I use curved grid to enhance accuracy

% b_grid grid for bond
% b_grid should be finer to see the difference in default probability under different productivity shocks at steady state
Nb = 20;
grid_b = curvspace(0,kMax,Nb,2)'; % To cover up as wide leverage level as possible

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
% laborFunction = @(a,k) ((k.^aalphaK * aalphaL).* a/W).^(1/(1-aalphaL));
% profitFunction = @(a,k,labor)a.* k.^aalphaK.* labor.^aalphaL - W*labor;
profitFunction = @(a,k)a.* k.^aalphaK.* ((k.^aalphaK * aalphaL).* a/W).^(aalphaL/(1-aalphaL)) - W * ((k.^aalphaK * aalphaL).* a/W).^(1/(1-aalphaL));
% nonDefaultFunction = @ (profit,k,bond)((1-ttaoC)*profit + (1-ddelta)*k > bond);
% isDefaultNextPeriod2DFunction =
% @(ia,mIsDefaultNextPeriod3D)(mIsDefaultNextPeriod3D(:,:,ia)); % never used below

% profit = a^(1/(1-aalphaL)) * (k.^aalphaK * (aalphaL * k.^aalphaK / W).^(aalphaL/(1-aalphaL)) - W * (aalphaL * k.^aalphaK/W).^(1/(1-aalphaL)));
% investmentFunction = @(k,kPrime,mIsDefaultNextPeriod)kPrime.*(1-mIsDefaultNextPeriod) - (1-ddelta)*k; %k_prime usually is k_grid
investmentFunction = @(k,kPrime)kPrime - (1-ddelta)*k; %k_prime usually is k_grid
% taxPaymentsFunction = @(k,bond,profit,RbMinus)ttaoC * (profit - ddelta*k - bond.* (RbMinus-1).* ((1-ttaoC)*profit + (1-ddelta)*k > bond)); % note the non-default indicator
taxPaymentsFunction = @(k,bond,profit,RbMinus)ttaoC * (profit - ddelta*k - bond * (RbMinus-1)); % note the non-default indicator
% taxPaymentsFunction = @(k,bond,profit,Rb)ttaoC * (profit - ddelta*k - bond * (Rb-1)); % without the non-default indicator

dividentFunction = @(profit,investment,bond,bondPrime,RbMinus,taxPayments,mIsDefaultNextPeriod)(profit - investment  ...
    + bondPrime.*(1-mIsDefaultNextPeriod) - RbMinus.* bond - taxPayments).*(1 + llambda * ((profit - investment  ...
    + bondPrime.*(1-mIsDefaultNextPeriod) - RbMinus.* bond - taxPayments) < 0)); % note the indicator function for issuance cost
% dividentFunction = @(profit,investment,bond,bondPrime,RbMinus,taxPayments)(profit - investment  ...
%     + bondPrime - RbMinus.* bond - taxPayments).*...
%     (1 + llambda * ((profit - investment + bondPrime - RbMinus.* bond - taxPayments) < 0)); % note the indicator function for issuance cost

%% a) Value Function Iteration when Rb = 1:01Rf
% Solve the Bellman equation for equity.
% Plot the optimal investment and default decision for the equity holders.

% State variables: k, b, z
% Control variables: k',b'

% Considerations: 
% aa) whether to default this period (using value function<0 this period to decide); 
% bb) if next period regardless of a', value function is all less than 0, then b' can't be borrowed this much,
% and thus the firm had to fund the investment by equity.

% Use multigrid method to speed up iteration
kGridLength           = [20]; % number of points in grid for capital
Nk = max(kGridLength);
kMin            = 0.000001;
kMax            = 7 * kSteadyState;
bMax            = 7 * kSteadyState;
grid_b = curvspace(0,bMax,Nb,2)';

% Required matrices and vectors
% Dimensionality is k,b,a,aMinus
    
kPolicyIndex = zeros(kGridLength(1),Nb,Na);
kPolicy = zeros(kGridLength(1),Nb,Na);
bPolicyIndex = zeros(kGridLength(1),Nb,Na);
bPolicy = zeros(kGridLength(1),Nb,Na);
% defaultPolicy = zeros(kGridLength(1),Nb,Na); % DON'T MAKE FOR DEFAULT A POLICY FUNCTION HERE! 
% Note default it not a % policy to make for the next period. default is a reaction in this period
% once productivity shock is realized and you see whether or not the value of the firm is larger than 0 or not.

mValue   = zeros(kGridLength(1),Nb,Na);
mValue0 = ones(kGridLength(1),Nb,Na); % initial guess


tic
for i=1:length(kGridLength)
    
    grid_k           = curvspace(kMin,kMax,kGridLength(i),2)';
    % Since the profitFunction takes so much time, let's calculate it all
    % at once to retrieve later
    mANkByNa = repmat(grid_a',kGridLength(i),1);
    mKNkByNa = repmat(grid_k,1,Na);
    profitNkByNa = profitFunction(mANkByNa,mKNkByNa); % Nk by Na
    
    tolerance = 0.00001;
    iteration = 0;
    distance = 100;

    kPrimeNkByNb = repmat(grid_k,1,Nb); % Nk*Nb matrix
    bondPrimeNkByNb = repmat(grid_b',kGridLength(i),1); % Nk*Nb matrix
    
    tic
    while distance > tolerance
        mIsDefaultNextPeriod = zeros(kGridLength(i),Nb); % k' b' entry is 1 if it will default for sure next period
        for ibPrime = 1:Nb
            for ikPrime = 1:kGridLength(i)
                mIsDefaultNextPeriod(ikPrime,ibPrime) = (sum((mValue0(ikPrime,ibPrime,:)>=0)) == 0);
            end
        end
        
%         tic
        for ia=1:Na
            a = grid_a(ia);

            for ik=1:kGridLength(i)
                k = grid_k(ik);        
                profit = profitNkByNa(ik,ia);% scalar
                investmentNkByNb = investmentFunction(k,kPrimeNkByNb); % Nk*Nb matrix

                for ib = 1:Nb
                    bond = grid_b(ib);
                    
                    if mValue0(ik,ib,ia) <0 % the firm will exit/default and choose no k' and b'
                        kPolicyIndex(ik,ib,ia) = 0;
                        kPolicy(ik,ib,ia) = 0;
                        bPolicyIndex(ik,ib,ia) = 0;
                        bPolicy(ik,ib,ia) = 0;
                    else % choose policy amongst k's b's
                        taxPayments = taxPaymentsFunction(k,bond,profit,RbConstant); % scalar
                        divident = dividentFunction(profit,investmentNkByNb,bond,bondPrimeNkByNb,RbConstant,taxPayments,mIsDefaultNextPeriod);% Nk*Nb matrix
%                         divident = (profit - investmentNkByNb + bondPrimeNkByNb？？？？？ - RbConstant.* bond - taxPayments).*...
%                             (1 + llambda * ((profit - investmentNkByNb + bondPrimeNkByNb - RbConstant.* bond - taxPayments) < 0));% Nk*Nb matrix

                        % Next we shall calculate the expected value
                        % tomorrow. Note we have to apply Max operator
                        % first and then take expectations - which means
                        % only when we reach that period will the firm make
                        % the default/nondefault decision.
                        mValueTomorrow = zeros(kGridLength(i),Nb,Na); % k',b',a'                            
                        for iaPrime = 1:Na % iterate over all possible states for tomorrow
                            mValueTomorrow(:,:,iaPrime) = m_a_prob(ia,iaPrime) * ...
                                (max(0,mValue0(:,:,iaPrime)).*(1-mIsDefaultNextPeriod)...
                                + max(0,repmat(mValue0(:,1,iaPrime),1,Nb)).*mIsDefaultNextPeriod);
                            
                        end
                        mExpectedValueTomorrow = sum(mValueTomorrow,3); % sum by the third dimension to get a Nk*Nb matrix

%                             mExpectedValueTomorrow = mExpectedValueTomorrowIfDefaultKnownForSureToday + mExpectedValueTomorrowIfDefaultKnownNotForSureToday;

                        x = divident + M * mExpectedValueTomorrow;

                        [rows,cols]=find(x==max(max(x)));

%                 if (1-ttaoC)*profit + (1-ddelta)*k <= bond
%                     value(ik,ib,ia) = 0;
%                 else
%                     value(ik,ib,ia) = x(max(rows),max(cols));
%                 end

                        kPolicyIndex(ik,ib,ia) = min(rows);
                        bPolicyIndex(ik,ib,ia) = min(cols);

                        kPolicy(ik,ib,ia) = grid_k(min(rows));
                        bPolicy(ik,ib,ia) = grid_b(min(cols));

                        mValue(ik,ib,ia) = max(max(x));

                    end
                end
            end
        end

%         toc
        %hold on
        %mesh(aa, kk, value')
        distance = sum(sum(sum(sum(abs(mValue(:,:,:,:)-mValue0(:,:,:,:))))));
        mValue0 = mValue;
        iteration = iteration + 1;

        if mod(iteration,5) == 0
            display("iteration =    " + iteration + "   difference =   " + distance )
        end
    end

    display("iteration =    " + iteration + "   difference =   " + distance + ". Converged")
    if i ~= length(kGridLength)
        mValue0 = interp1(grid_k, mValue,linspace(kMin, kMax, kGridLength(i+1))');% 这里不知道linspace后面要不要加'变为列向量。试试吧
        mValue  = mValue0;
%         kPolicy = interp1(grid_k,kPolicy,linspace(kMin, kMax, kGridLength(i+1)));
        kPolicy         = zeros(kGridLength(i+1),Nb,Na);
        kPolicyIndex    = zeros(kGridLength(i+1),Nb,Na);
        bPolicy         = zeros(kGridLength(i+1),Nb,Na);
        bPolicyIndex    = zeros(kGridLength(i+1),Nb,Na);
    end
    
end


toc
% save('valuePrevious','mValue')
save('resultA','mValue','kPolicy','bPolicy','kPolicyIndex','bPolicyIndex')

figure(1);
[bb,kk]=meshgrid(grid_b, grid_k);
mesh(bb, kk, mValue(:,:,1));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, mValue(:,:,ia));    
end

title('Value Under Different Shocks','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])

zlabel('Value','interpreter','latex')
savefig('q3a_value_3D')

figure(2)
mesh(bb, kk, kPolicy(:,:,1));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, kPolicy(:,:,ia));    
end

title('Policy $k^\prime$ Under Different Shocks','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])

zlabel('$k^\prime$','interpreter','latex')
savefig('q3a_kPolicy_3D')


figure(3)
mesh(bb, kk, bPolicy(:,:,1));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, bPolicy(:,:,ia));    
end

title('Policy $b^\prime$ Under Different Shocks','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])

zlabel('$bond^\prime$','interpreter','latex')
savefig('q3a_bPolicy_3D')

% Plot the optimal investment and default decision for the equity holders.
mK3D=repmat(grid_k,1,Nb,Na);% Nk*Nb*Na matrix
investmentPolicy = kPolicy - (1-ddelta) * mK3D;

figure(4)
mesh(bb, kk, investmentPolicy(:,:,1));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, investmentPolicy(:,:,ia));    
end

title('Policy $investment^\prime$ Under Different Shocks','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('$investment^\prime$','interpreter','latex')
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])

savefig('q3a_investmentPolicy_3D')

figure(5)
mIsDefaultToday = (mValue<0);
mesh(bb, kk, mIsDefaultToday(:,:,1));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, mIsDefaultToday(:,:,ia));    
end

title('Whether to Default Today Under Different Shocks','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('default','interpreter','latex')
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])

savefig('q3a_IsDefaultToday_3D')

% Compute the stationary distribution to see how good is the range of the
% grid

distributionStationary0 = (1/(Nk*Nb*Na))*ones(Nk,Nb,Na);
distance=100;
tolerance=0.0000000001;
iteration=0;

while distance>tolerance
    distributionStationary1 = zeros(Nk,Nb,Na);
    for ia=1:Na
%         for iaMinus=1:Na
        for ib=1:Nb
            for ik=1:Nk
                ikPrime = max(1,kPolicyIndex(ik,ib,ia));          
                ibPrime = max(1,bPolicyIndex(ik,ib,ia)); 

                prob = distributionStationary0(ik,ib,ia);
                for iaPrime=1:Na
                    prob_aPrime = prob*m_a_prob(ia,iaPrime);
                    distributionStationary1(ikPrime,ibPrime,iaPrime) = distributionStationary1(ikPrime,ibPrime,iaPrime) + prob_aPrime;
                end
                    
                
            end
        end
    end
    
    distance=sum(sum(sum(abs(distributionStationary0-distributionStationary1))));
    distributionStationary0 = distributionStationary1;
    iteration = iteration + 1;
end

% Plot the distribution

figure(6);
[bb,kk]=meshgrid(grid_b, grid_k);
mesh(bb, kk, distributionStationary0(:,:,1));
for ia = 2:Na
    hold on
    mesh(bb, kk, distributionStationary0(:,:,ia));
end

title('Stationary Distribution','interpreter','latex');
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('Probability Mass','interpreter','latex')
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])

savefig('q3a_stationary_distribution_3D')



%% b) Probability of Default Next Period
% b)  Compute and plot the probability of default next period, 
% conditional on the value of the shocks today p(z; b'; k')
mDefaultProbability = zeros(Nk,Nb,Na);% k',b',a
for ibPrime = 1:Nb
    for ikPrime = 1:Nk
        for ia = 1:Na
            mDefaultProbability(ikPrime,ibPrime,ia) = sum(reshape((mValue(ikPrime,ibPrime,:)<0),1,Na).*m_a_prob(ia,:));
        end
    end
end

figure(7)
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
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])

savefig('q3b_default_prob_3D')

%% Required Rate of Return of Bonds
% c) Now use the conditional default probability to compute the required rate of return
% by bondholders, Rb(b'; k'; z) that ensures they make 0 profits. Assume for simplicity
% bondholders get paid 0 upon default.

mRf = Rf*ones(Nk,Nb,Na);
mCarryOnProbability = 1-mDefaultProbability;
mIsDefaultNextPeriod3D = (mDefaultProbability==1); % (b'; k'; z)
% mRb = min(Rf/mCarryOnProbability,1000000000);
mRb = Rf/mCarryOnProbability;
% mRb(~isfinite(mRb))=NaN;% replace the inf with NaN

%% Value Function Iteration given Rb = mRb(k',b',z).
% d) Solve the Bellman equation for the equity holders taking as given this new function
% for the required rate of return by bondholders, Rb(·). Plot this new value function
% against the one found in section a).

% Now we have to keep track of a-.
% State variables: k,b,a,a-
% Control variables: k',b'

% Use multigrid method to speed up iteration
% kGridLength           = [20]; % number of points in grid for capital
% Nk = max(kGridLength);
% kMin            = 0.000001;
% kMax            = 5 * kSteadyState;
% bMax            = 5 * kSteadyState;

% grid_b = curvspace(0,bMax,Nb,2)';

% Required matrices and vectors
% Dimensionality is k,b,a,aMinus
    
kPolicyIndex = zeros(kGridLength(1),Nb,Na,Na);
kPolicy = zeros(kGridLength(1),Nb,Na,Na);
bPolicyIndex = zeros(kGridLength(1),Nb,Na,Na);
bPolicy = zeros(kGridLength(1),Nb,Na,Na);
% defaultPolicy = zeros(kGridLength(1),Nb,Na); % DON'T MAKE FOR DEFAULT A POLICY FUNCTION HERE! 
% Note default it not a % policy to make for the next period. default is a reaction in this period
% once productivity shock is realized and you see whether or not the value of the firm is larger than 0 or not.

mValue   = zeros(kGridLength(1),Nb,Na,Na);
mValue0 = ones(kGridLength(1),Nb,Na,Na); % initial guess


tic
for i=1:length(kGridLength)
    
    grid_k           = curvspace(kMin,kMax,kGridLength(i),2)';
    % Since the profitFunction takes so much time, let's calculate it all
    % at once to retrieve later
    mANkByNa = repmat(grid_a',kGridLength(i),1);
    mKNkByNa = repmat(grid_k,1,Na);
    profitNkByNa = profitFunction(mANkByNa,mKNkByNa); % Nk by Na
    
    tolerance = 0.00001;
    iteration = 0;
    distance = 100;

    kPrimeNkByNb = repmat(grid_k,1,Nb); % Nk*Nb matrix
    bondPrimeNkByNb = repmat(grid_b',kGridLength(i),1); % Nk*Nb matrix
    
    tic
    while distance > tolerance
        % Now I don't have to calculte mIsDefaultNextPeriod since I already compute that.
        
%         tic
        for ia=1:Na
            a = grid_a(ia);
            mIsDefaultNextPeriod2D = mIsDefaultNextPeriod3D(:,:,ia); % Nk*Nb matrix

            for ik=1:kGridLength(i)
                k = grid_k(ik);        
                profit = profitNkByNa(ik,ia);% scalar
                investmentNkByNb = investmentFunction(k,kPrimeNkByNb); % Nk*Nb matrix

                for ib = 1:Nb
                    bond = grid_b(ib);
                    
                    for iaMinus = 1:Na
                        RbMinus = mRb(ik,ib,iaMinus); % scalar
%                         aMinus = grid_a(iaMinus);
                    
                        if mValue0(ik,ib,ia,iaMinus) <0 % the firm will exit/default and choose no k' and b'
                            kPolicyIndex(ik,ib,ia,iaMinus) = 0;
                            kPolicy(ik,ib,ia,iaMinus) = 0;
                            bPolicyIndex(ik,ib,ia,iaMinus) = 0;
                            bPolicy(ik,ib,ia,iaMinus) = 0;
                        else % choose policy amongst k's b's
                            taxPayments = taxPaymentsFunction(k,bond,profit,RbMinus); % scalar
                            divident = dividentFunction(profit,investmentNkByNb,bond,bondPrimeNkByNb,RbMinus,taxPayments,mIsDefaultNextPeriod2D);% Nk*Nb matrix
    %                         divident = (profit - investmentNkByNb + bondPrimeNkByNb？？？？？ - RbConstant.* bond - taxPayments).*...
    %                             (1 + llambda * ((profit - investmentNkByNb + bondPrimeNkByNb - RbConstant.* bond - taxPayments) < 0));% Nk*Nb matrix

                            % Next we shall calculate the expected value
                            % tomorrow. Note we have to apply Max operator
                            % first and then take expectations - which means
                            % only when we reach that period will the firm make
                            % the default/nondefault decision.
                            mValueTomorrow = zeros(kGridLength(i),Nb,Na); % k',b',a'                            
                            for iaPrime = 1:Na % iterate over all possible states for tomorrow
                                mValueTomorrow(:,:,iaPrime) = m_a_prob(ia,iaPrime) * ...
                                    (max(0,mValue0(:,:,iaPrime)).*(1-mIsDefaultNextPeriod2D)...
                                    + max(0,repmat(mValue0(:,1,iaPrime),1,Nb)).*mIsDefaultNextPeriod2D);
                            
                            end
                            mExpectedValueTomorrow = sum(mValueTomorrow,3); % sum by the third dimension to get a Nk*Nb matrix

                            x = divident + M * mExpectedValueTomorrow;

                            [rows,cols]=find(x==max(max(x)));

    %                 if (1-ttaoC)*profit + (1-ddelta)*k <= bond
    %                     value(ik,ib,ia) = 0;
    %                 else
    %                     value(ik,ib,ia) = x(max(rows),max(cols));
    %                 end

                            kPolicyIndex(ik,ib,ia,iaMinus) = min(rows);
                            bPolicyIndex(ik,ib,ia,iaMinus) = min(cols);

                            kPolicy(ik,ib,ia,iaMinus) = grid_k(min(rows));
                            bPolicy(ik,ib,ia,iaMinus) = grid_b(min(cols));

                            mValue(ik,ib,ia,iaMinus) = max(max(x));

                        end
                    end
                end
            end
        end

%         toc
        %hold on
        %mesh(aa, kk, value')
        distance = sum(sum(sum(sum(abs(mValue(:,:,:,:)-mValue0(:,:,:,:))))));
        mValue0 = mValue;
        iteration = iteration + 1;

        if mod(iteration,5) == 0
            display("iteration =    " + iteration + "   difference =   " + distance )
        end
    end

    display("iteration =    " + iteration + "   difference =   " + distance + ". Converged")
    if i ~= length(kGridLength)
        mValue0 = interp1(grid_k, mValue,linspace(kMin, kMax, kGridLength(i+1))');% 这里不知道linspace后面要不要加'变为列向量。试试吧
        mValue  = mValue0;
%         kPolicy = interp1(grid_k,kPolicy,linspace(kMin, kMax, kGridLength(i+1)));
        kPolicy         = zeros(kGridLength(i+1),Nb,Na,Na);
        kPolicyIndex    = zeros(kGridLength(i+1),Nb,Na,Na);
        bPolicy         = zeros(kGridLength(i+1),Nb,Na,Na);
        bPolicyIndex    = zeros(kGridLength(i+1),Nb,Na,Na);
    end
    
end


toc
% save('valuePrevious','mValue')
save('resultD','mValue','kPolicy','bPolicy','kPolicyIndex','bPolicyIndex')


figure(8);
[bb,kk]=meshgrid(grid_b, grid_k);
mesh(bb, kk, mValue(:,:,1,round((Na+1)/2)));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, mValue(:,:,ia,round((Na+1)/2)));    
end

title('Value Under Different Shocks given mean $z^{-}$','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('Value','interpreter','latex')
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])

zlim([-max(max(max(max(mValue)))),max(max(max(max(mValue))))])
% legend('low productivity','median productivity','high productivity') % no use in mesh
savefig('q3d_value_3D')

figure(9)
mesh(bb, kk, kPolicy(:,:,1,round((Na+1)/2)));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, kPolicy(:,:,ia,round((Na+1)/2)));    
end

title('Policy $k^\prime$ Under Different Shocks given mean $z^{-}$','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('$k^\prime$','interpreter','latex')
% legend('low productivity','median productivity','high productivity') % no use in mesh
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])

savefig('q3d_kPolicy_3D')


figure(10)
mesh(bb, kk, bPolicy(:,:,1,round((Na+1)/2)));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, bPolicy(:,:,ia,round((Na+1)/2)));    
end

title('Policy $b^\prime$ Under Different Shocks given mean $z^{-}$','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('$bond^\prime$','interpreter','latex')
% legend('low productivity','median productivity','high productivity') % no use in mesh
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])
savefig('q3d_bPolicy_3D')

% Compute the stationary distribution to see how good is the range of the
% grid
distributionStationary0 = (1/(Nk*Nb*Na*Na))*ones(Nk,Nb,Na,Na);
distance=100;
tolerance=0.0000000001;
iteration=0;

while distance>tolerance
    distributionStationary1 = zeros(Nk,Nb,Na,Na);
    for ia=1:Na
        for iaMinus=1:Na
            for ib=1:Nb
                for ik=1:Nk
                    ikPrime = max(1,kPolicyIndex(ik,ib,ia,iaMinus));          
                    ibPrime = max(1,bPolicyIndex(ik,ib,ia,iaMinus)); 

                    prob = distributionStationary0(ik,ib,ia,iaMinus);
                    for iaPrime=1:Na
                        prob_aPrime = prob*m_a_prob(ia,iaPrime);
                        distributionStationary1(ikPrime,ibPrime,iaPrime,ia) = distributionStationary1(ikPrime,ibPrime,iaPrime,ia) + prob_aPrime;
                    end
                    
                end
            end
        end
    end
    
    distance=sum(sum(sum(sum(abs(distributionStationary0-distributionStationary1)))));
    distributionStationary0 = distributionStationary1;
    iteration = iteration + 1;
end

% Plot the distribution

figure(11);
[bb,kk]=meshgrid(grid_b, grid_k);
mesh(bb, kk, distributionStationary0(:,:,1,1));
for iaMinus = 2:Na
    for ia = 2:Na
        hold on
        mesh(bb, kk, distributionStationary0(:,:,ia,iaMinus));
    end
end
title('Stationary Distribution','interpreter','latex');
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('Probability Mass','interpreter','latex')
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])

savefig('q3d_stationary_distribution_3D')

% open q3d_value_3D.fig
% open q3a_value_3D.fig

% Plot this new value function against the one found in section a)

figure(12);
[bb,kk]=meshgrid(grid_b, grid_k);

subplot(1,2,1)
load('resultA.mat')
mesh(bb, kk, mValue(:,:,1));

for ia = 2:Na
    hold on;
    mesh(bb, kk, mValue(:,:,ia));    
end
title('Value Functions a)','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])
zlim([-max(max(max(max(mValue)))),max(max(max(max(mValue))))])
zlabel('Value','interpreter','latex')

subplot(1,2,2)
load('resultD.mat')
mesh(bb, kk, mValue(:,:,1));
for ia = 2:Na
    hold on;
    mesh(bb, kk, mValue(:,:,ia));    
end
title('Value Functions d)','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
ylim([min(grid_k),max(grid_k)])
xlim([min(grid_b),max(grid_b)])
zlim([-max(max(max(max(mValue)))),max(max(max(max(mValue))))])
zlabel('Value','interpreter','latex')
savefig('q3d_value_3D_resultA_and_resultD')

% As we can see from the picture, the difference is only at the region where 
% capital is extremely low and bond is extremely high.

% This makes sense, because in a model where interest rate on risky bond
% does not depend on the default risk but is exogenously given, there will
% hardly be negative values because you won't pay very high interest rate
% even though your default probability is almost 1, or exactly one. That
% case, if you don't default, you just pay a relatively low interest rate 
% and if you default because the value would be negative if you didn't, 
% the supposed negative value cannot be too negative because the interest 
% rate is just 1.01 times the riskfree rate.

% On the other hand, in a model where interest rate on risky bond does
% depend on the default risk, once you ask for a very risky debt, the
% interest rate is so high that once there is a low productivity shock next
% period, you are more willing to default because if you chose not to, the
% value would be very very negative because you would have to pay a super
% high interest rate for the high risk debt you asked for last period.