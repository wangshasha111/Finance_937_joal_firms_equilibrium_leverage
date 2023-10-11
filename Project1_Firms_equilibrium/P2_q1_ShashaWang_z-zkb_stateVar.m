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
kMin = 0.00001;

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
mIsDefaultNextPeriod = (mDefaultProbability==1);

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

%% Functions
% create functions for convenience
% laborFunction = @(a,k) ((k.^aalphaK * aalphaL).* a/W).^(1/(1-aalphaL));
% profitFunction = @(a,k,labor)a.* k.^aalphaK.* labor.^aalphaL - W*labor;
profitFunction = @(a,k)a.* k.^aalphaK.* ((k.^aalphaK * aalphaL).* a/W).^(aalphaL/(1-aalphaL)) - W * ((k.^aalphaK * aalphaL).* a/W).^(1/(1-aalphaL));
nonDefaultFunction = @ (profit,k,bond)((1-ttaoC)*profit + (1-ddelta)*k > bond);
% isDefaultNextPeriod2DFunction =
% @(ia,mIsDefaultNextPeriod3D)(mIsDefaultNextPeriod3D(:,:,ia)); % never used below

% profit = a^(1/(1-aalphaL)) * (k.^aalphaK * (aalphaL * k.^aalphaK / W).^(aalphaL/(1-aalphaL)) - W * (aalphaL * k.^aalphaK/W).^(1/(1-aalphaL)));
% investmentFunction = @(k,kPrime,mIsDefaultNextPeriod)kPrime.*(1-mIsDefaultNextPeriod) - (1-ddelta)*k; %k_prime usually is k_grid
investmentFunction = @(k,kPrime)kPrime - (1-ddelta)*k; %k_prime usually is k_grid
taxPaymentsFunction = @(k,bond,profit,RbMinus)ttaoC * (profit - ddelta*k - bond.* (RbMinus-1).* ((1-ttaoC)*profit + (1-ddelta)*k > bond)); % note the non-default indicator
% taxPaymentsFunction = @(k,bond,profit,Rb)ttaoC * (profit - ddelta*k - bond * (Rb-1)); % note the non-default indicator
% taxPaymentsFunction = @(k,bond,profit,Rb)ttaoC * (profit - ddelta*k - bond * (Rb-1)); % without the non-default indicator

dividentFunction = @(profit,investment,bond,bondPrime,RbMinus,taxPayments,mIsDefaultNextPeriod)(profit - investment  ...
    + bondPrime.*(1-mIsDefaultNextPeriod) - RbMinus.* bond - taxPayments).*(1 + llambda * ((profit - investment  ...
    + bondPrime.*(1-mIsDefaultNextPeriod) - RbMinus.* bond - taxPayments) < 0)); % note the indicator function for issuance cost

%% c) Value Function Iteration
% Solve the Bellman equation for the equity holders taking as given the
% function for the required rate of return by bondholders, Rb

% Use multigrid method to speed up iteration
kGridLength           = [15]; % number of points in grid for capital
kMin            = 0.000001;
kMax            = 10 * kSteadyState;
grid_b = curvspace(0,kMax,Nb,2)';

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
    
    grid_k           = curvspace(kMin,kMax,kGridLength(i),2)';
    % Since the profitFunction takes so much time, let's calculate it all
    % at once to retrieve later
    mANkByNa = repmat(grid_a',kGridLength(i),1);
    mKNkByNa = repmat(grid_k,1,Na);
    profitNkByNa = profitFunction(mANkByNa,mKNkByNa); % Nk by Na

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
    mIsDefaultNextPeriod = (mDefaultProbability==1);
    mCarryOnProbability = 1-mDefaultProbability;
    mRb = min(Rf/mCarryOnProbability,1000000);
    
    tolerance = 0.00001;
    iteration = 0;
    distance = 100;

    kPrime = repmat(grid_k,1,Nb); % Nk*Nb matrix
    bondPrime = repmat(grid_b',kGridLength(i),1); % Nk*Nb matrix
    
    % transition matrix for (k',b',a, k)->k'_actually
%     kPrimeIfDefaultKnownForSureToday = zeros(kGridLength(i),Nb,kGridLength(i),Na); % k',b',k
%     for ik = 1:kGridLength(i)
%         comparison = abs(grid_k-grid_k(ik)*(1-ddelta));
%         ind_kPrime = find(comparison == min(comparison));
%         for ia = 1:Na
%         kPrimeIfDefaultKnownForSureToday(:,:,ik) = grid_k(ind_kPrime).*;
%         
%     for ibPrime = 1:Nb
%         for ikPrime = 1:kGridLength(i)
%             comparison = abs(grid_k-grid_k(ik)*(1-ddelta)); % See how much capital there will actually be left tomorrow
%             ind_kPrime = find(comparison == min(comparison));
%             valueTomorrowIfDefaultKnownForSureToday(ikPrime,ibPrime,:) = isDefaultNextPeriod(ikPrime,ibPrime) * reshape(m_a_prob(ia,:),1,1,3).* reshape(value0(ind_kPrime,1,:,ia),1,1,3);
%         end
%     end
%     mExpectedValueTomorrowIfDefaultKnownForSureToday = sum(valueTomorrowIfDefaultKnownForSureToday,3); % sum by the third dimension to get a Nk*Nb matrix


    tic
    while distance > tolerance
%         tic
        for ia=1:Na
            a = grid_a(ia);
            isDefaultNextPeriod = mIsDefaultNextPeriod(:,:,ia); % Nk*Nb matrix

            for ik=1:kGridLength(i)
                k = grid_k(ik);        
%                 labor = laborFunction(a,k);
%                 profit = profitFunction(a,k,labor);
%                 profitFunction(a,k); % scalar
                profit = profitNkByNa(ik,ia);% scalar
                investment = investmentFunction(k,kPrime); % Nk*Nb matrix

                for iaMinus = 1:Na
                    RbMinus = mRb(:,:,iaMinus); % Nk*Nb matrix
                    aMinus = grid_a(iaMinus);

                    for ib = 1:Nb
                        bond = grid_b(ib);

                        if (1-ttaoC)*profit + (1-ddelta)*k <= bond % if default this period
                            value (ik,ib,ia,iaMinus)=0; % You stop operating the firm and stop choosing next period k' and b'
                            % nothing will happen to policy function index matrix or function matrix - entry remains 0
                            
%                             kPolicyIndex(ik,ib,ia,iaMinus) = 1;
%                             bPolicyIndex(ik,ib,ia,iaMinus) = 1;
% 
%                             kPolicy(ik,ib,ia,iaMinus) = grid_k(1);
%                             bPolicy(ik,ib,ia,iaMinus) = grid_b(1);

                        else % if not default this period

                            taxPayments = taxPaymentsFunction(k,bond,profit,RbMinus); % Nk*Nb matrix
                            divident = dividentFunction(profit,investment,bond,bondPrime,RbMinus,taxPayments,isDefaultNextPeriod); % Nk*Nb matrix
                             
%                             valueTomorrowIfDefaultKnownForSureToday = zeros(kGridLength(i),Nb,Na); % k',b',a'
%                             for iaPrime = 1:Na
%                                 aPrime = grid_a(iaPrime);
%                                 kPrime = %Nk*Nb matrix
% %                                 laborPrime = laborFunction(aPrime,kPrime);
% %                                 profitPrime = profitFunction(aPrime,kPrime,laborPrime);
% %                                 profitPrime = profitFunction(aPrime,kPrime); % Nk*Nb
%                                 profitPrime = repmat(profitNkByNa(:,iaPrime),1,Nb);% Nk*Nb %也需要retrieve
%                                 
%                             valueTomorrowIfDefaultKnownForSureToday(:,:,iaPrime) = value0(:,:,iaPrime,ia).*isDefaultNextPeriod.* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime > bondPrime);
%                             for ibPrime = 1:Nb
%                                 for ikPrime = 1:kGridLength(i)
%                                     comparison = abs(grid_k-grid_k(ik)*(1-ddelta)); % See how much capital there will actually be left tomorrow
%                                     ind_kPrime = find(comparison == min(comparison));
%                                     valueTomorrowIfDefaultKnownForSureToday(ikPrime,ibPrime,:) = isDefaultNextPeriod(ikPrime,ibPrime) * reshape(m_a_prob(ia,:),1,1,3).* reshape(value0(ind_kPrime,1,:,ia),1,1,3);
%                                 end
%                             end
%                             mExpectedValueTomorrowIfDefaultKnownForSureToday = sum(valueTomorrowIfDefaultKnownForSureToday,3); % sum by the third dimension to get a Nk*Nb matrix
                          
                            
                            valueTomorrow = zeros(kGridLength(i),Nb,Na); % k',b',a'                            
                            for iaPrime = 1:Na % iterate over all possible states for tomorrow
                                aPrime = grid_a(iaPrime);
                                profitPrime = repmat(profitNkByNa(:,iaPrime),1,Nb);% Nk*Nb
                                valueTomorrow(:,:,iaPrime) = m_a_prob(ia,iaPrime) * ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime > bondPrime)...% 需要考虑default之后value为0
                                    .*( value0(:,:,iaPrime,ia).*(1-isDefaultNextPeriod)...% if known today for sure that NO default next period
                                    + repmat(value0(:,1,iaPrime,ia),1,Nb).*isDefaultNextPeriod); % if known today for sure WILL default next period
                            end
                            mExpectedValueTomorrow = sum(valueTomorrow,3); % sum by the third dimension to get a Nk*Nb matrix
                            
%                             mExpectedValueTomorrow = mExpectedValueTomorrowIfDefaultKnownForSureToday + mExpectedValueTomorrowIfDefaultKnownNotForSureToday;
                            
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
save('resultC','value','kPolicy','bPolicy','kPolicyIndex','bPolicyIndex')

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

% Note now the setting has changed a little bit from question c). 
% All we have to change is the value function continuation in the bellman
% equation, where we don't set value of default to zero. We set if to value
% when k and b are 0.

% Notice this is the only difference from c). So I basically copied the
% codes for c) and made the change accordingly. 
% 1) This period, if the firm % default, instead of setting the value to 0,
% I set it to value at k=0.0000001 and b=0, and policy function entry is 
% set to equal to the entry when k=0.0000001 and b=0.
% 2) Next Period, if the firm % default, instead of setting the valueTomorrow
% to 0, I set it to value at k=0.0000001 and b=0.
% I use two kinds of endogenous grids to solve the problem - (a-,a,k,b) and
% (a,n). For the first method, we can easily make both the first and second
% adjustment, but for the second method where we only have today's networth
% information, we can only make the second adjustment and assume that
% today's firms are nondefaulters and the 0 networth denotes defaulters.

% Use multigrid method to speed up iteration
kGridLength           = [15]; % number of points in grid for capital
Nk = max(kGridLength);
kMin            = 0.000001;
kMax            = 10 * kSteadyState;
grid_b = curvspace(0,kMax,Nb,2)';

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
    
    grid_k           = curvspace(kMin,kMax,kGridLength(i),2)';
    
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
    mIsDefaultNextPeriod = (mDefaultProbability==1);
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
            isDefaultNextPeriod = mIsDefaultNextPeriod(:,:,ia); % Nk*Nb matrix
            
            for ik=1:kGridLength(i)
                k = grid_k(ik);        
%                 labor = laborFunction(a,k);
%                 profit = profitFunction(a,k,labor);        
%                 profit = profitFunction(a,k); 
                profit = profitNkByNa(ik,ia);
                investment = investmentFunction(k,kPrime); % Nk*Nb matrix

                for iaMinus = 1:Na
                    RbMinus = mRb(:,:,iaMinus); % Nk*Nb matrix
                    aMinus = grid_a(iaMinus);

                    for ib = 1:Nb
                        bond = grid_b(ib);

                        if (1-ttaoC)*profit + (1-ddelta)*k <= bond % if default this period
                            value (ik,ib,ia,iaMinus)=value0(1,1,ia,iaMinus); % after restructuring, you continue to run the firm at square 1 with 0 capital and 0 bond
                            
                            kPolicyIndex(ik,ib,ia,iaMinus) = kPolicyIndex(1,1,ia,iaMinus) ; % policy function is the same as when you are at square 1
                            bPolicyIndex(ik,ib,ia,iaMinus) = bPolicyIndex(1,1,ia,iaMinus) ;

                            kPolicy(ik,ib,ia,iaMinus) = grid_k(kPolicyIndex(1,1,ia,iaMinus));
                            bPolicy(ik,ib,ia,iaMinus) = grid_b(bPolicyIndex(1,1,ia,iaMinus));

                        else % if not default this period

                            taxPayments = taxPaymentsFunction(k,bond,profit,RbMinus); % Nk*Nb matrix
                            divident = dividentFunction(profit,investment,bond,bondPrime,RbMinus,taxPayments,isDefaultNextPeriod); % Nk*Nb matrix

                            valueTomorrow = zeros(kGridLength(i),Nb,Na);% k',b',a'
                            
                            for iaPrime = 1:Na % iterate over all possible states for tomorrow
                                aPrime = grid_a(iaPrime);
                                profitPrime = repmat(profitNkByNa(:,iaPrime),1,Nb);% Nk*Nb
                               
%                                 valueTomorrow(:,:,iaPrime) =  m_a_prob(ia,iaPrime) * ... % 需要考虑default之后value不是为0，而是为set bond and capital to 0的value
%                                     (value0(:,:,iaPrime,ia).* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime > bondPrime)... % if not default next period
%                                     + value0(1,1,iaPrime,ia).* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime <= bondPrime));% if default next period
%                                 
                                valueTomorrow(:,:,iaPrime) =  m_a_prob(ia,iaPrime) *( ...% 需要考虑default之后value不是为0，而是为set bond and capital to 0的value
                                    (1-isDefaultNextPeriod).*(...% if known today not for sure that default will happen next period
                                    value0(:,:,iaPrime,ia).* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime > bondPrime)... 
                                    + value0(1,1,iaPrime,ia).* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime <= bondPrime))...
                                    +isDefaultNextPeriod.* ... % if known today for sure WILL default next period
                                    (repmat(value0(:,1,iaPrime,ia),1,Nb).* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime > bondPrime)... % next period if not default 
                                    + value0(1,1,iaPrime,ia).* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime <= bondPrime)));% next period if default 
                            end
                            mExpectedValueTomorrow = sum(valueTomorrow,3); % sum by the third dimension to get a Nk*Nb matrix

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
save('resultD','value','kPolicy','bPolicy','kPolicyIndex','bPolicyIndex')

figure(6);
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
savefig('q1d_value_3D')

figure(7)
mesh(bb, kk, kPolicy(:,:,1,round((Na+1)/2)));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, kPolicy(:,:,ia,round((Na+1)/2)));    
end

title('Policy $k^\prime$ Under Different Shocks given mean $z^{-}$','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('$k^\prime$','interpreter','latex')
savefig('q1d_kPolicy_3D')


figure(8)
mesh(bb, kk, bPolicy(:,:,1,round((Na+1)/2)));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, bPolicy(:,:,ia,round((Na+1)/2)));    
end

title('Policy $b^\prime$ Under Different Shocks given mean $z^{-}$','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('$bond^\prime$','interpreter','latex')
savefig('q1d_bPolicy_3D')


%% Stationary Distribution
% Compute the stationary distribution of firms. 

distributionStationary0 = (1/(Nk*Nb*Na*Na))*ones(Nk,Nb,Na,Na);
distance=100;
tolerance=0.000001;
iteration=0;

while distance>tolerance
    distributionStationary1 = zeros(Nk,Nb,Na,Na);
    for ia=1:Na
        for iaMinus=1:Na
            for ib=1:Nb
                for ik=1:Nk
                    ikPrime = kPolicyIndex(ik,ib,ia,iaMinus);          
                    ibPrime = bPolicyIndex(ik,ib,ia,iaMinus);

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
figure(9);
[bb,kk]=meshgrid(grid_b, grid_k);
aMinusDescription = ["low","medium","high"];

for iaMinus = 1:Na
    subplot(1,Na,iaMinus);
    mesh(bb, kk, distributionStationary0(:,:,1,iaMinus));
    
    for ia = 2:Na
        hold on;
        mesh(bb,kk,distributionStationary0(:,:,ia,iaMinus));
    end
    title(['Distribution $z^{-}$ ',aMinusDescription(iaMinus)],'interpreter','latex');
    ylabel('Capital Stock $k^\prime$','interpreter','latex')
    xlabel('Debt $b^\prime$','interpreter','latex')
    zlabel('Probability Mass','interpreter','latex')
end
    
savefig('q1d_stationary_distribution_3D')

%% Table
% Use the invariant distribution to construct
% a table reporting the cross-sectional average values of:
% (1) probability of default, p(·);
% (2) required return on risky bonds, Rb(·);
% (3) leverage ratio, b/k;
% (4) investment to capital ratio, i/k.
% (5) fraction of firms issuing equity;

% To do this question, I created 4D array to reduce the layer of loops.

mK4D=repmat(grid_k,1,Nb,Na,Na);% Nk*Nb*Na*Na matrix
mBond4D=repmat(grid_b',Nk,1,Na,Na); % Nk*Nb*Na*Na matrix

mA4D=repmat(grid_a',Na,1,Nb,Nk);
mA4D=permute(mA4D,[4,3,2,1]);% transform/reallocate the dimension to get a Nk*Nb*Na*Na matrix

mAMinus4D=repmat(grid_a,1,Na,Nb,Nk);
mAMinus4D=permute(mAMinus4D,[4,3,2,1]);% transform/reallocate the dimension to get a Nk*Nb*Na*Na matrix

%% (1) probability of default, p(·);

% mLabor4D = laborFunction(mA4D,mK4D); % Nk*Nb matrix
% mProfit4D = profitFunction(mA4D,mK4D,mLabor4D);% Nk*Nb matrix
mProfit4D = profitFunction(mA4D,mK4D);% Nk*Nb matrix

mDefault4D = 1-nonDefaultFunction(mProfit4D,mK4D,mBond4D);% Nk*Nb 0-1 matrix

defaultProbability = sum(sum(sum(sum(mDefault4D.*distributionStationary0))));
fprintf('Average Default Probability is %2.10f\n', defaultProbability);

% mK = repmat(grid_k,1,Nb); % Nk*Nb matrix
% mBond = repmat(grid_b',Nk,1); % Nk*Nb matrix

% defaultAll = zeros(Nk,Nb,Na,Na);
% for iaMinus = 1:Na
%     for ia = 1:Na
%         a = grid_a(ia);
%         mLabor = laborFunction(a,mK); % Nk*Nb matrix
%         mProfit = profitFunction(a,mK,mLabor);% Nk*Nb matrix
%         mDefault = 1-nonDefaultFunction(mProfit,mK,mBond);% Nk*Nb 0-1 matrix
%         defaultAll(:,:,ia,iaMinus) = mDefault;
%     end
% end
% 
% defaultProbability = sum(sum(sum(sum(defaultAll.*distributionStationary0))));


%% (2) required return on risky bonds, Rb(·);
% riskyBondReturn = sum(sum(sum(sum((min(1000000000,Rf/(1-mDefault4D))).*distributionStationary0))));
riskyBondReturn = sum(sum(sum(sum((min(10000000000000,Rf/(1-mDefault4D))).*((mDefault4D~=1).*distributionStationary0)))));
fprintf('Required return on risky bonds on average is %2.8f\n', riskyBondReturn);
% As we can see, corporate bonds are basically riskfree.

% riskyBondReturn = Rf/(1-defaultProbability);
% fprintf('Required return on risky bonds is %2.8f\n', riskyBondReturn);

%% (3) leverage ratio, b/k;
leverageRatio = sum(sum(sum(sum(mBond4D ./ mK4D   .*distributionStationary0))));
fprintf('Average leverage ratio is %2.8f\n', leverageRatio);

% leverageRatioAll = zeros(Nk,Nb,Na,Na);
% for iaMinus = 1:Na
%     for ia = 1:Na
%         a = grid_a(ia);
%         leverageRatioAll(:,:,ia,iaMinus) = mBond./mK;
%     end
% end
% 
% leverageRatio = sum(sum(sum(sum(leverageRatioAll.*distributionStationary0))));
% fprintf('Average leverage ratio is %2.8f\n', leverageRatio);

%% (4) investment to capital ratio, i/k.
investment2Capital=sum(sum(sum(sum((kPolicy./mK4D + 1 - ddelta).*distributionStationary0))));
fprintf('Average investment to capital ratio is %2.8f\n', investment2Capital);

%% (5) fraction of firms issuing equity;
    % mIsDefaultNextPeriod4D = repmat(mIsDefaultNextPeriod,Na,1);行不通
mIsDefaultNextPeriod4D = zeros(kGridLength(1),Nb,Na,Na);
for iaMinus=1:Na
    aMinus = grid_a(iaMinus);
    mIsDefaultNextPeriod4D(:,:,:,iaMinus)=mIsDefaultNextPeriod;
end

mRbMinus4D = zeros(kGridLength(1),Nb,Na,Na);
for ia = 1:Na
    a = grid_a(ia);
    for iaMinus = 1:Na
        aMinus = grid_a(iaMinus);
        mRbMinus4D(:,:,ia,iaMinus)=mRb(:,:,iaMinus);
    end
end
    
mInvestment4D = investmentFunction(mK4D,mK4D);
mTaxPayments4D = taxPaymentsFunction(mK4D,mBond4D,mProfit4D,mRbMinus4D);
mDivident4D = dividentFunction(mProfit4D,mInvestment4D,mBond4D,mBond4D,mRbMinus4D,mTaxPayments4D,mIsDefaultNextPeriod4D);

mIsIssuingEquity = (mDivident4D<0);
fractionOfFirmsIssuingEquity = sum(sum(sum(sum(  mIsIssuingEquity .*distributionStationary0))));
fprintf('The fraction of firms issuing equity is %2.8f\n', fractionOfFirmsIssuingEquity);

T = table(defaultProbability,riskyBondReturn,leverageRatio,investment2Capital,fractionOfFirmsIssuingEquity)
%     defaultProbability    riskyBondReturn       leverageRatio      investment2Capital    fractionOfFirmsIssuingEquity
%     __________________    ________________    _________________    __________________    ____________________________
% 
%   3.56994771855542e-09    1.01010100649501    0.539685179800697     2.00309413226414          0.030296293479109      

% Nk=20, Nb=20, Na=5
%     defaultProbability    riskyBondReturn    leverageRatio    investment2Capital    fractionOfFirmsIssuingEquity
%     __________________    _______________    _____________    __________________    ____________________________
% 
%         3.5699e-09            1.0101            0.65464             1.989                     0.030296      

