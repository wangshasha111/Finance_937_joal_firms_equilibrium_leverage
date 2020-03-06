%%
tic
clc
clear all
close all

%% Parameters
M = 0.95;
delta = 0.1;
alpha_k = 0.3;
alpha_l = 0.65;
W = 2;
theta = 0.5;
rho = 0.95;
sigma = 0.02;

k_stead = 1;
r = 1/M-1;
alpha = alpha_k/(1-alpha_l);
a_bar = k_stead^(1-alpha)*(r+delta)/alpha;
z_bar = a_bar^(1-alpha_l)/(1-alpha_l)^(1-alpha_l)/(alpha_l/W)^alpha_l;
logzbar = log(z_bar);
% logzbar = 0.5;

Ent = 0.025;

%% Variables

std_coeff_z_plus = [0.25:0.25:1.5,2,2.5,3,4];
std_coeff_z_minus = - std_coeff_z_plus(end:-1:1);
std_coeff_z = [std_coeff_z_minus,0,std_coeff_z_plus];

Nkg = 501;
precision=10^-8;

Nzg = length(std_coeff_z);

 
%% Making Grids for ``z"  and ``k"


logz_vec = logzbar + std_coeff_z*sigma/(1-rho^2)^0.5 ;
Pz=zeros(Nzg,Nzg);

logz_vec_exted = [-Inf,logz_vec,Inf];
cdf_trunc = @(x,mean,std) max(min((cdf('Normal',x,mean,std) - cdf('Normal',mean-4*std,mean,std))/(cdf('Normal',mean+4*std,mean,std) - cdf('Normal',mean-4*std,mean,std)),1),0);   
% cdf_trunc = @(x,mean,std) (cdf('Normal',x,mean,std) - cdf('Normal',mean-4*std,mean,std))/(cdf('Normal',mean+4*std,mean,std) - cdf('Normal',mean-4*std,mean,std));   

for i=1:Nzg
    
 Pz(i,:) = cdf_trunc((logz_vec_exted(2:end-1)+logz_vec_exted(3:end))/2,rho*logz_vec(i)+(1-rho)*logzbar,sigma) - cdf_trunc((logz_vec_exted(1:end-2)+logz_vec_exted(2:end-1))/2,rho*logz_vec(i)+(1-rho)*logzbar,sigma);
    
end

z_vec = exp(logz_vec);


k_vec_max = 25;
k_vec_min = 0.01;


k_vec = exp(linspace(log(k_vec_min),log(k_vec_max),Nkg));
 
 %% Calibration
error_stat = 1;
max_stat_eval = 1200;

f_range = [0,0.05];

while error_stat>5e-4
f = mean(f_range);



%% Instant. Payoff Functions
 pai = @(z,k) k.^(alpha_k/(1-alpha_l)) .* (alpha_l*z/W).^(1/(1-alpha_l))*(1/alpha_l-1)*W - f;
 phi = @(i,k)  theta*(i./k-delta).^2.*k;

 

%% Individual Problem Solver
k_mat = repmat(reshape(k_vec,[Nkg,1,1]),[1,Nzg,Nkg]);
z_mat = repmat(reshape(z_vec,[1,Nzg,1]),[Nkg,1,Nkg]);
kp_mat = repmat(reshape(k_vec,[1,1,Nkg]),[Nkg,Nzg,1]);
i_mat = kp_mat - (1-delta) * k_mat;
div =  pai(z_mat,k_mat) - i_mat - phi(i_mat,k_mat);


maxeval=1000;

initguess_V = zeros(Nkg,Nzg);
V=initguess_V;
for l=1:maxeval 
    
   V_prev=V;
     
      EV=V*Pz';     
      V_temp = div + M * repmat(reshape(max(EV',0),[1,Nzg,Nkg]),[Nkg 1 1]); 
      
      Policy_stay = EV>= 0;
      
      [V,pos]=max(V_temp,[],3);
        
        
     diff=max(max(abs(V-V_prev)));
    if diff<precision
        disp('Policy function Convergence achieved')
        break
    end
    if l==maxeval
        disp('MaxEval achieved, no convergence')
    end
%     if mod(l,10)==0 || l==1
%         disp('Iteration   Error Ratio :')
%         disp([l,diff/precision])
%     end
    
end


%% Policy Functions
k_mat = repmat(reshape(k_vec,[Nkg,1]),[1,Nzg]);
z_mat = repmat(reshape(z_vec,[1,Nzg]),[Nkg,1]);

Policy_k = k_mat(pos);
Policy_i = Policy_k - (1-delta)*k_mat;
Policy_d =  pai(z_mat,k_mat) - Policy_i - phi(Policy_i,k_mat);



%% Stationary Distribution

% E=1;
% Policy_stay_vec = Policy_stay(:);
% pos_vec = pos(:);
% temp = repmat(0:Nzg-1,Nkg,1);
% pos_vec = pos_vec + temp(:)*Nkg;
% T = kron (Pz' , eye(Nkg));
% 
% temp = eye(Nkg*Nzg);
% temp = temp(1:Nkg*Nzg,pos_vec);
% temp = repmat(eye(Nkg),Nzg,Nzg)*temp;
% T = T*temp;
% T = T./(repmat(sum(T),Nkg*Nzg,1));
% 
% [station_dist,station_eig_dist] = eig(T);
% station_dist = station_dist(:,1);
% 
% if abs(station_eig_dist(1)-1)>1e-5 
%     'Error'
% end


[station_state_z,station_eig_z] = eig(Pz');
station_state_z = station_state_z(:,1);
station_state_z = station_state_z/sum(station_state_z);

if abs(station_eig_z(1)-1)>1e-5 
    'Error'
end



total_mass = zeros(1,max_stat_eval);
mu = zeros(Nkg,Nzg)/Nkg/Nzg ;

disp('current f is:')
disp(f)
        disp('total mass evolution:') 

for i = 1:max_stat_eval
    
    total_mass(i) = sum(sum(mu));
    if mod(i,200)==0
        disp(total_mass(i))
    end
    
    new_mu = zeros(Nkg,Nzg);
    for j = 1:Nzg
        for k=1:Nkg
          new_mu(pos(k,j),:) = new_mu(pos(k,j),:) + (Policy_stay(k,j)*mu(k,j))*Pz(j,:);       
        end
    end
    
    new_mu(1,:) = new_mu(1,:) + Ent * station_state_z';
    
         diff=max(max(abs(new_mu/sum(sum(new_mu))-mu/sum(sum(mu)))));
         
%     if diff<precision/Nkg/Nzg/10
%         disp('Convergence achieved')
%         break
%     end
%     if l==maxeval
%         disp('MaxEval achieved, no convergence in distribution')
%     end
%     if mod(i,10)==0 || i==1
%         disp('Iteration On Distribution Error Ratio :')
%         disp([i,diff/precision*Nkg*Nzg])
%     end
    
    Exit_rate = sum(sum(mu.*(1-Policy_stay)))/sum(sum(mu));

%     disp('exit rate:')
%     disp(Exit_rate)
    
    
    mu = new_mu;


end




%% adjust f
disp('derived exit_rate is : ')
disp(Exit_rate)


error_stat = abs(Exit_rate-Ent);

if Exit_rate<Ent
    f_range(1) = f;
    disp('too low, increase f !!')
end
if Exit_rate>Ent
    f_range(2) = f; 
    disp('too high, decrease f !!')
end




end



%%




Entry_Cost = sum(station_state_z'.*V(1,:));
% 
disp('% of Exit Rate : ')
disp(Exit_rate*100)
% 
% 
save('data_Q1.mat')











toc