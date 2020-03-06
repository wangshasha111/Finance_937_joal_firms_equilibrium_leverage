%%

tic

clc
close all


try
    rho = rho;
    lambda = lambda;
catch
    clear 
   rho = 0.7; 
   lambda = 0.025; 
end

 



%% Parameters
M = 0.99;
delta = 0.1;
alpha_k = 0.3;
alpha_l = 0.6;
W = 2;

sigma = 0.05;
logzbar = 1.25; % is set arbitrary -- to have readable scales.

alpha = alpha_k/(1-alpha_l);
A = @(z) (alpha_l*z/W).^(1/(1-alpha_l))*(1/alpha_l-1)*W;
Z  = @(a) ( a/W/(1/alpha_l-1) ).^(1-alpha_l)*W/alpha_l;



tau_c = 0.15;

R_f = 1/M;



Nzg = 5;
Nkg = 40;
Nbg = 20;

precision = 1e-8;
maxeval = 10000;

%% Making Grids for ``z"  and ``k"
logz_vec = logzbar - (Nzg-1)/2*sigma/(1-rho^2)^0.5 + (0:Nzg-1)*sigma/(1-rho^2)^0.5;
Pz=zeros(Nzg,Nzg);

logz_vec_exted = [-Inf,logz_vec,Inf];
cdf_trunc = @(x,mean,std) max(min((cdf('Normal',x,mean,std) - cdf('Normal',mean-3*std,mean,std))/(cdf('Normal',mean+3*std,mean,std) - cdf('Normal',mean-3*std,mean,std)),1),0);   


for i=1:Nzg
    
 Pz(i,:) = cdf_trunc((logz_vec_exted(2:end-1)+logz_vec_exted(3:end))/2,rho*logz_vec(i)+(1-rho)*logzbar,sigma) - cdf_trunc((logz_vec_exted(1:end-2)+logz_vec_exted(2:end-1))/2,rho*logz_vec(i)+(1-rho)*logzbar,sigma);
    
end

z_vec = exp(logz_vec);
a_vec = A(z_vec);


k_vec_max = ( a_vec((Nzg+1)/2) / delta )^(1/(1-alpha));
%k_vec_min = ( A_vec(1) / (delta+R_f) )^(1/(1-alpha));
k_vec_min = 0.01;
% k_s = ( alpha*a_vec((Nzg+1)/2) / (delta+R_f) )^(1/(1-alpha)); % this is static optimization level of k which is about 1 due to choice of logzbar
% k_vec = unique([linspace(k_vec_min,k_s,Nkg/2+1),exp(linspace(log(k_s),log(k_vec_max),Nkg/2))]);
k_vec = exp(linspace(log(k_vec_min),log(k_vec_max),Nkg));


b_vec_max = k_vec_max*1.5;
% b_s = k_s/3*2 ;
b_vec_min = k_vec_min;
% b_vec = unique([linspace(b_vec_min,b_s,Nbg/2+1),exp(linspace(log(b_s),log(b_vec_max),Nbg/2))]);
b_vec = exp(linspace(log(b_vec_min),log(b_vec_max),Nbg));

Num_state = Nzg^2;
Num_choice = Nkg*Nbg;


%% Q. a) Function Default Probability , Plot

pai = @(A,k) k.^alpha .* A; 

%%% k' and b' are row vector inputs of same size here. A_index is a vector
%%% of integers poiting to elements in vector grid of productivity.
P_Def = @(A_index,bp,kp) (Pz(A_index,:)*(repmat(a_vec',1,length(bp)).*repmat(kp.^alpha,length(a_vec),1)*(1-tau_c) + repmat(kp,length(a_vec),1)*(1-delta) < repmat(bp,length(a_vec),1)))';
% P_Def = @(A_index,bp,kp) (normcdf( (repmat(log(Z(max(bp - (1-delta)*kp,0)./(kp.^alpha)/(1-tau_c))),length(A_index),1) - ( (1-rho)*logzbar + rho * repmat(logz_vec(A_index)',1,length(bp)) ) )/sigma  ))';

 
plot_counter = 1;
for now_kp_plot = [0.2,2,8]
    
    bp_range_plot = linspace(0,15,1000);
    kp_range_plot = now_kp_plot*ones(size(bp_range_plot));
    
figure(plot_counter)
hold on
plot(repmat(bp_range_plot,3,1)',P_Def([1,(Nzg+1)/2,Nzg],bp_range_plot,kp_range_plot),'linewidth',1.5)
title('Part I: a) \it Default Probability','interpreter','Latex','fontsize',15)
xlabel('$b''$','interpreter','Latex','fontsize',14)
ylabel('$P^{Def.}(z,b'',k'')$','interpreter','Latex','fontsize',14)
legend({[repmat(['$k'' = ',num2str(now_kp_plot),', ~z='],3,1),['low ';'med.';'high'],repmat('$',3,1)]},'interpreter','Latex','location','best','fontsize',13)
set(gca,'YTick',0:0.2:1)
box on
plot_counter = plot_counter+1;

end



%% Q. b) breakeven rate for investors:


P_b = @(A_index,bp,kp) (1-P_Def(A_index,bp,kp))/R_f;
R_e_b = @(A_index,bp,kp) min(1./P_b(A_index,bp,kp),1e20);

% importnat note: P_D = 0 implies R_e_b = +inf. there might be then issues in
% working with this. if your'e code is sensitive to this you can fix this by adding an infiniteimal number to
% repay probability or take a min with some arbitrary large number, larger than maximum taxshield benefit.




%% Q. c) Belman Equation, given price for debt


% first define flow of payoff; inputs should be comparable in dimension!!!
div = @ (prev_A_index,A,k,b,k_p,b_p) (1-tau_c)*pai(A,k) - k_p +(1-delta*(1-tau_c))*k + b_p - tau_c*b -(1-tau_c)*b.*repmat(kron(R_e_b(prev_A_index(1,1:Nzg:end,1),b(:,1,1)',k(:,1,1)'),ones(1,Nzg)),[1,1,Num_choice]);


% I am going to build 3d arrays
% put cur. k and b (as exog states for current period) in first element
% put prev. A and current A (as exog states) in second elemnt
% put k' and b' (as choice vars.) in third element.
% this is really inefficient; it's better to put evey exogenous vars. in
% one dimension, so that we could convert matrix to sparse. more kron.
% products are needed to do that and the code gets really unreadable.


prev_A_index_mat = repmat(reshape(kron(1:Nzg,ones(1,Nzg)),[1,Num_state,1]),[Num_choice,1,Num_choice]); % note: prev A index takes values of 1 to Nzg
A_mat = repmat(reshape(kron(ones(1,Nzg),a_vec),[1,Num_state,1]),[Num_choice,1,Num_choice]);

k_mat = repmat(reshape(kron(ones(1,Nbg),k_vec),[Num_choice,1,1]),[1,Num_state,Num_choice]);
b_mat = repmat(reshape(kron(b_vec,ones(1,Nkg)),[Num_choice,1,1]),[1,Num_state,Num_choice]);

kp_mat = repmat(reshape(kron(ones(1,Nbg),k_vec),[1,1,Num_choice]),[Num_choice,Num_state,1]);
bp_mat = repmat(reshape(kron(b_vec,ones(1,Nkg)),[1,1,Num_choice]),[Num_choice,Num_state,1]);

% transition probability matrix for exog states: prev A index and A 
Pz_mat = repmat(repmat(Pz,1,Nzg).*kron(eye(Nzg),ones(1,Nzg)),Nzg,1);

% now construct flow payoff 3d array
F = div(prev_A_index_mat,A_mat,k_mat,b_mat,kp_mat,bp_mat);

issu_indx = F<0;
F = F.*(1 + lambda.*issu_indx );

Default = pai(A_mat,k_mat)*(1-tau_c) + k_mat * (1-delta) < b_mat;


initguess_V = zeros(Num_choice,Num_state); % not a good guess! but let's see weather we can converge even with crazy initial guess (robust enough...)
% initguess_V = 1/(1-M) * F(:,:,Nkg/2); % not a good guess! but...
V=initguess_V;
for l=1:maxeval 
    
   V_prev=V;
     
      EV=V*Pz_mat'; 
%        EV(isnan(EV)) = 0;  % if never going to a state with -inf value  (prob = 0) then assign zero value to that event on expectation
           
     V_temp = F + M * repmat(reshape(EV',[1,Num_state,Num_choice]),[Num_choice 1 1]);
     % if default, liquidate firm asset, no continuation value.
     V_temp(Default & (V_temp>=-1e16)) = 0 ; % limited liability if happened to default right now... but if everyone knew the firm would default by sure (R_e = 1e20) would'nt be allowed to borrow at all in previous period
     
      if sum(sum(sum(isnan(V_temp))))>0
          disp('warning: V has nan values!!')
      end
     
%         F_quit = F(:,:,1);
%       V(~stay_index(:,:,1))= F_quit(~stay_index(:,:,1)) ; 
%         V(isinf(V))=-1e16;
      [V,pos]=max(V_temp,[],3);



        error = V-V_prev;
%         error(isnan(error)) = 0;
     diff=max(max(abs(error)));
    if diff<precision
        disp('Convergence achieved')
        break
    end
    if l==maxeval
        disp('MaxEval achieved, no convergence')
    end
    if mod(l,10)==0 || l==1
        disp('Iteration ,  Error Ratio :')
        disp(l)
        disp(diff/precision)
    end
    
end



% plot to see how value function looks like.

figure(plot_counter)
plot(1:Nkg*Nbg,V(:,Nzg*(Nzg-1)/2+1:(Nzg-1)/2:Nzg*(Nzg+1)/2),'linewidth',1.5)
title(['\it Part 1, Q. c):  Value Function'],'interpreter','Latex','fontsize',15)
xlabel('$[kb]$','interpreter','Latex','fontsize',14)
ylabel('$V(z,[kb])$','interpreter','Latex','fontsize',14)
legend({[repmat('prev. z  = med. \& curr. $z=',3,1),['min';'med';'max'],repmat('$',3,1)]},'interpreter','Latex','location','best','fontsize',13)
box on
plot_counter = plot_counter + 1;


% Policy Functions

stay_index = ~Default;
% there is essentially no policy to stay, it is ad-hoc-ly pinned down by
% current stat!
Policy_stay = stay_index(:,:,1);

Policy_k = k_mat(pos);
Policy_k(~Policy_stay) = nan;

Policy_b = b_mat(pos);
Policy_b(~Policy_stay) = nan;


Policy_i = Policy_k - (1-delta)*k_mat(:,:,1);
Policy_i(~Policy_stay) = nan;


Policy_issue = nan(size(V));
for i=1:Num_choice 
    for j=1:Num_state
Policy_issue(i,j) = issu_indx(i,j,pos(i,j));
    end
end
Policy_issue(~Policy_stay) = nan;


% plot to see how policy functions look like.
%%% note: policy function is not assigned if the firm quits. iterestinlgy
%%% the firm will quit if current z realization is low. will issue maximum
%%% possible debt on grid if curr. z is median to enjoy tax shield, even
%%% though it knows it will default and quit by sure next period, but will
%%% not do so if current z realization is high since it wants to stay and
%%% enjoy high profit flow of operation.

figure(plot_counter)
plot(1:Nkg*Nbg,Policy_b(:,Nzg*(Nzg-1)/2+1:(Nzg-1)/2:Nzg*(Nzg+1)/2),'linewidth',1.5)
title(['\it Part 1, Q. c):  Policy Function'],'interpreter','Latex','fontsize',15)
xlabel('$[kb]$','interpreter','Latex','fontsize',14)
ylabel('$b''(z,[kb])$','interpreter','Latex','fontsize',14)
legend({[repmat('prev. z  = med. \& curr. $z=',3,1),['min';'med';'max'],repmat('$',3,1)]},'interpreter','Latex','location','best','fontsize',13)
box on
plot_counter = plot_counter + 1;

figure(plot_counter)
plot(1:Nkg*Nbg,Policy_k(:,Nzg*(Nzg-1)/2+1:(Nzg-1)/2:Nzg*(Nzg+1)/2),'linewidth',1.5)
title(['\it Part 1, Q. c):  Policy Function'],'interpreter','Latex','fontsize',15)
xlabel('$[kb]$','interpreter','Latex','fontsize',14)
ylabel('$k''(z,[kb])$','interpreter','Latex','fontsize',14)
legend({[repmat('prev. z  = med. \& curr. $z=',3,1),['min';'med';'max'],repmat('$',3,1)]},'interpreter','Latex','location','best','fontsize',13)
box on
plot_counter = plot_counter + 1;


%% Q. d) Stationary Dist
% [station_state_z,station_eig_z] = eig(Pz');
% station_state_z = station_state_z(:,1);
% station_state_z = station_state_z/sum(station_state_z);

mu = ones(Num_choice,Num_state)/Num_choice/Num_state ; % this is crazy initial guess, but let's see weather algorithm is sensitive to initial guess or not.

for i = 1:maxeval
    
    new_mu = zeros(Num_choice,Num_state);
    for j = 1:Num_state
        for k=1:Num_choice
          new_mu(pos(k,j),:) = new_mu(pos(k,j),:) + mu(k,j)*Pz_mat(j,:)*Policy_stay(k,j);
          new_mu(1,:) = new_mu(1,:)+ mu(k,j)*Pz_mat(j,:)*(1-Policy_stay(k,j));
        end
    end
        
         diff=max(max(abs(new_mu-mu)));
         
    if diff<precision/Num_choice/Num_state/1e8
        disp('Convergence achieved')
        break
    end
    if l==maxeval
        disp('MaxEval achieved, no convergence in distribution')
    end
    if mod(i,10)==0 || i==1
        disp('Iteration On Distribution , Error Ratio :')
        disp(i)
        disp(diff/precision*1e8*Num_choice*Num_state)
    end
    
    
    mu = new_mu;
    
end



% plot the stationary distibution

figure(plot_counter)
surf(1:Nzg,1:Nkg*Nbg,mu(:,Nzg*(Nzg-1)/2+1:1:Nzg*(Nzg+1)/2),'linewidth',1.5)
title('\it Part 1, Q. d):  Stationary Distribution','interpreter','Latex','fontsize',15)
xlabel('$[z]$','interpreter','Latex','fontsize',14)
ylabel('$[kb]$','interpreter','Latex','fontsize',14)
zlabel('$\mu([kb],z)$','interpreter','Latex','fontsize',14)
legend({'prev. [z] = median'},'interpreter','Latex','location','best','fontsize',13)
box on
plot_counter = plot_counter + 1;

figure(plot_counter)
surf(1:Nzg,1:Nkg*Nbg,mu(:,Nzg*(Nzg-1)+1:1:Nzg*Nzg),'linewidth',1.5)
title('\it Part 1, Q. d):  Stationary Distribution','interpreter','Latex','fontsize',15)
xlabel('$[z]$','interpreter','Latex','fontsize',14)
ylabel('$[kb]$','interpreter','Latex','fontsize',14)
zlabel('$\mu([kb],z)$','interpreter','Latex','fontsize',14)
legend({'prev. [z] = max'},'interpreter','Latex','location','best','fontsize',13)
box on
plot_counter = plot_counter + 1;

figure(plot_counter)
surf(1:Nzg,1:Nkg*Nbg,mu(:,1:1:Nzg),'linewidth',1.5)
title('\it Part 1, Q. d):  Stationary Distribution','interpreter','Latex','fontsize',15)
xlabel('$[z]$','interpreter','Latex','fontsize',14)
ylabel('$[kb]$','interpreter','Latex','fontsize',14)
zlabel('$\mu([kb],z)$','interpreter','Latex','fontsize',14)
legend({'prev. [z] = min'},'interpreter','Latex','location','best','fontsize',13)
box on
plot_counter = plot_counter + 1;







%% Outputs

disp('---------------------')
disp('------ statistics -------------')

b_mat = b_mat(:,:,1);
k_mat = k_mat(:,:,1);
A_mat = A_mat(:,:,1);
pai = pai(A_mat,k_mat);

mu_current = mu;
mu_current = mu_current/sum(sum(mu_current));

mu_stay = mu;
mu_stay(~Policy_stay) = 0;
ave_Prob_Default = 1-sum(sum(mu_stay))



R_b_mat = kron(R_e_b(prev_A_index_mat(1,1:Nzg:end,1),b_mat(:,1)',k_mat(:,1)'),ones(1,Nzg));
ave_R_D_temp = R_b_mat.*mu_current;
ave_R_D_temp(isnan(ave_R_D_temp))=0;
ave_R_e_D = sum(sum(ave_R_D_temp))


ave_Leverage_temp = b_mat./k_mat.*mu_current;
ave_Leverage = sum(sum(ave_Leverage_temp))

ave_iok_temp = Policy_i./k_mat.*mu_current;
ave_iok = nansum(nansum(ave_iok_temp))

ave_Issue = nansum(nansum(Policy_issue.*mu_current))


if min(min(V.*mu))<-1e16
    'Error: Positive math firm with -Inf Value'
else
    'Sounds good: No active firm with -Inf Values'
    V(V<-1e16)=0;
end


%% Saving Data

% save('Data_P1.mat','A_mat','k_mat','b_mat','pai','V','mu','Policy_stay','Policy_issue','Policy_i','Policy_k','Policy_b');





toc


