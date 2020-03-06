
%%

tic



clc
clear all
close all



%% Parameters
M = 0.99;
delta = 0.1;
alpha_k = 0.3;
alpha_l = 0.6;
W = 2;

rho = 0.7;
sigma = 0.05;
logzbar = 1.25;

alpha = alpha_k/(1-alpha_l);
A = @(z) (alpha_l*z/W).^(1/(1-alpha_l))*(1/alpha_l-1)*W;


lambda = 0.025;
tau_c = 0.15;

R_f = 1/M;


Nzg = 7;
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
A_vec = A(z_vec);


k_vec_max = ( A_vec((Nzg+1)/2) / delta )^(1/(1-alpha));
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

Num_state = Nzg;
Num_choice = Nkg*Nbg;


% %% 3D Vectorized
% 
% A_mat = repmat(reshape(A_vec,[1,Num_state,1]),[Num_choice,1,Num_choice]);
% 
% k_mat = repmat(reshape(kron(ones(1,Nbg),k_vec),[Num_choice,1,1]),[1,Num_state,Num_choice]);
% b_mat = repmat(reshape(kron(b_vec,ones(1,Nkg)),[Num_choice,1,1]),[1,Num_state,Num_choice]);
% 
% kp_mat = repmat(reshape(kron(ones(1,Nbg),k_vec),[1,1,Num_choice]),[Num_choice,Num_state,1]);
% bp_mat = repmat(reshape(kron(b_vec,ones(1,Nkg)),[1,1,Num_choice]),[Num_choice,Num_state,1]);
% 
% Pz_mat = Pz;

%% (full) 2D Vectorized

A_mat = kron(ones(1,Num_choice),A_vec);
A_mat_full = sparse(repmat(reshape(A_mat,[1,Num_choice*Num_state]),[Num_choice,1]));

k_mat_prim = kron(ones(1,Nbg),k_vec);
b_mat_prim = kron(b_vec,ones(1,Nkg));
k_mat = kron(k_mat_prim,ones(1,Num_state));
b_mat = kron(b_mat_prim,ones(1,Num_state));
k_mat_full = sparse(repmat(reshape(k_mat,[1,Num_choice*Num_state]),[Num_choice,1]));
b_mat_full = sparse(repmat(reshape(b_mat,[1,Num_choice*Num_state]),[Num_choice,1]));

kp_mat_prim = kron(ones(1,Nbg),k_vec);
bp_mat_prim = kron(b_vec,ones(1,Nkg));
kp_mat_full = sparse(repmat(reshape(kp_mat_prim,[Num_choice,1]),[1,Num_choice*Num_state]));
bp_mat_full = sparse(repmat(reshape(bp_mat_prim,[Num_choice,1]),[1,Num_choice*Num_state]));

temp = eye(Num_choice);
temp = sparse(repmat(temp(:),[1,Num_choice]));
Pz_mat_prime = kron(  temp , Pz');

%%

 
R_b =  1.01 * R_f * ones(Num_choice,Num_state*Num_choice);
R_e_b = R_b*(1-tau_c) + tau_c;
P_b = 1./R_e_b;

num_big_Iter = 2;
for big_Iter = 1:num_big_Iter


    

pai = @(A,k) k.^alpha .* A; 
div = @ (A,k,b,k_p,b_p) (1-tau_c)*pai(A,k) - k_p +(1-delta*(1-tau_c))*k + P_b.*b_p - b;




%% Belman Equation


F = div(A_mat_full,k_mat_full,b_mat_full,kp_mat_full,bp_mat_full);

issu_indx = F<0;
F = sparse(F.*(1 + lambda.*issu_indx ));


V = zeros(1,Num_choice*Num_state);
for l=1:maxeval 
    
   V_prev=V;
     
   EV = kron(sparse(eye(Num_choice)),V) * Pz_mat_prime; % trsutt me[?] this works :)
      V_temp = F + M * EV;
%       V_temp(V_temp<bp_mat_full)=-inf; %firms cannot borrow money more than its value, since it would then default by sure...
        
      [V,pos]=max(V_temp,[],1);
          V = full(V);
     diff=max(abs(V-V_prev));
    if diff<precision
%         disp('Convergence achieved')
        break
    end
    if l==maxeval
        disp('MaxEval achieved, no convergence')
    end
    if mod(l,10)==0 || l==1
        disp('Iteration  , Error Ratio :')
        disp(l)
        disp(diff/precision)
    end
    
end

Default_indx = V<0;

% first method to recalibrate price of debt:
P_e_Def = kron(sparse(eye(Num_choice)),Default_indx) * Pz_mat_prime;  % debt holder take the next-period-default-behavior of previous loop step as given 
R_b = R_f./(1-P_e_Def);
R_e_b = R_b*(1-tau_c) + tau_c;
P_b = 1./R_e_b;
% P_b_new = nan(Num_choice,Num_state*Num_choice);
% for pricing_indx = 1:Num_choice*Num_state
%     P_b_new(:,pricing_indx) = P_b; % price of debt is independent of initial state b,k
% end
% P_b = P_b_new;


%%% second method to recalibrate price of debt, based on return rate on
%%% debt
% P_e_Def = zeros(Num_choice,Num_state);
% for jjjj=1:Num_state
%     Pos_JJ = zeros(1,Num_choice^2);
%     Pos_JJ(pos(:,jjjj)'+(0:Num_choice-1)*Num_choice)=1;
%     Pos_JJ = reshape(Pos_JJ,[Num_choice,Num_choice])';
%     P_e_Def(:,jjjj) = Pos_JJ * Default_indx * (Pz_mat(jjjj,:))';
% end
% P_b =  (1-P_e_Def)/R_f;
% R_e_b = 1./P_b';
% R_b = R_e_b(1:Nzg,:);
% R_b = kron(R_b,ones(Nzg,1));
% R_b = R_b';
% R_b = repmat(R_b,[1,1,Num_choice]); % debtholds are pricing the debt based on the previous loop step optimal oplicy of equity holders, which is a function of initial state. note: b'(k,b,z) k'(k,b,z). equity holders may rather deviate actually and choose another level of debt.




% disp('First Round Finished.')
% disp('Is there anyone who expected to not default?')
% disp(any(P_e_Def(:)<1))
% 
% disp('Is there anyone who expected to Default?')
% disp(any(P_e_Def(:)>0))

if big_Iter > 1
%     disp(max(abs(R_b_prev(:)-R_b(:))))
    if max(abs(R_e_b_prev(:)-R_e_b(:)))<1e-4*R_f
        'YEEESSSS'
        break
    end
end

R_e_b_prev = R_e_b;


end








toc


