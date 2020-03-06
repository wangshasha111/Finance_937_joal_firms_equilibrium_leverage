
rho = 0.7;
lambda = 0.025;

run('Part_I.m')


%% Q.a)

tbl_1 = table(Policy_i(:)./k_mat(:),V(:)./k_mat(:),pai(:)./k_mat(:),'VariableNames',{'i_over_k','Q','pi_over_k'});
lm_1 = fitlm(tbl_1,'i_over_k~Q+pi_over_k','Weights',mu_stay(:))


tbl_2 = table(b_mat(:)./k_mat(:),V(:)./k_mat(:),pai(:)./k_mat(:),log(k_mat(:)),'VariableNames',{'b_over_k','Q','pi_over_k','log_k'});
lm_2 = fitlm(tbl_2,'b_over_k ~ Q + pi_over_k + log_k','Weights',mu_stay(:))





%% Q.b)



index_iterr = 1;


rho_range = [0.5,0.6,0.7,0.8,0.9];
lambda_range = [0,0.0125,0.025,0.05,0.1];

for i_estim_indx=1:length(rho_range)
    for j_estim_indx=1:length(lambda_range)
       rho = rho_range(i_estim_indx);
       lambda = lambda_range(j_estim_indx);
       
        run('Part_I.m')

tbl_1 = table(Policy_i(:)./k_mat(:),V(:)./k_mat(:),pai(:)./k_mat(:),'VariableNames',{'i_over_k','Q','pi_over_k'});
lm_1b_{i_estim_indx,j_estim_indx} = fitlm(tbl_1,'i_over_k~Q+pi_over_k','Weights',mu_stay(:));

tbl_2 = table(b_mat(:)./k_mat(:),V(:)./k_mat(:),pai(:)./k_mat(:),log(k_mat(:)),'VariableNames',{'b_over_k','Q','pi_over_k','log_k'});
lm_2b_{i_estim_indx,j_estim_indx} = fitlm(tbl_2,'b_over_k ~ Q + pi_over_k + log_k','Weights',mu_stay(:));

        index_iterr = index_iterr+1;
        
    end
end



beta_1=zeros(5,5);
beta_2=zeros(5,5);
gamma_1=zeros(5,5);
gamma_2=zeros(5,5);
gamma_3=zeros(5,5);

beta_1_std=zeros(5,5);
beta_2_std=zeros(5,5);
gamma_1_std=zeros(5,5);
gamma_2_std=zeros(5,5);
gamma_3_std=zeros(5,5);

for i_estim_indx=1:length(rho_range)
    for j_estim_indx=1:length(lambda_range)
        
        disp('rho = ')
        disp(rho_range(i_estim_indx))  
        
        disp('lambda = ')
        disp(lambda_range(j_estim_indx))
      
        
        beta_1(i_estim_indx,j_estim_indx) = lm_1b_{i_estim_indx,j_estim_indx}.Coefficients{2,1};
        beta_2(i_estim_indx,j_estim_indx) = lm_1b_{i_estim_indx,j_estim_indx}.Coefficients{3,1};

        beta_1_std(i_estim_indx,j_estim_indx) = lm_1b_{i_estim_indx,j_estim_indx}.Coefficients{2,2};
        beta_2_std(i_estim_indx,j_estim_indx) = lm_1b_{i_estim_indx,j_estim_indx}.Coefficients{3,2};
        
        gamma_1(i_estim_indx,j_estim_indx) = lm_2b_{i_estim_indx,j_estim_indx}.Coefficients{2,1};
        gamma_2(i_estim_indx,j_estim_indx) = lm_2b_{i_estim_indx,j_estim_indx}.Coefficients{3,1};
        gamma_3(i_estim_indx,j_estim_indx) = lm_2b_{i_estim_indx,j_estim_indx}.Coefficients{4,1};

        gamma_1_std(i_estim_indx,j_estim_indx) = lm_2b_{i_estim_indx,j_estim_indx}.Coefficients{2,2};
        gamma_2_std(i_estim_indx,j_estim_indx) = lm_2b_{i_estim_indx,j_estim_indx}.Coefficients{3,2};
        gamma_3_std(i_estim_indx,j_estim_indx) = lm_2b_{i_estim_indx,j_estim_indx}.Coefficients{4,2};
        
        
%         disp('beta_1 , beta_2 , gamma_1 , gamma_2 , gamma_3 = ')
%         
%         disp(beta_1(i,j))
%         disp(beta_2(i,j))
%         disp(gamma_1(i,j))
%         disp(gamma_2(i,j))
%         disp(gamma_3(i,j))
                
        
    end
end



X_target = [0.004,0.2,-0.15,-0.4,0.05];


error = zeros(5,5);
for i_estim_indx=1:length(rho_range)
    for j_estim_indx=1:length(lambda_range)
        
        X_temp = [beta_1(i_estim_indx,j_estim_indx),beta_2(i_estim_indx,j_estim_indx),gamma_1(i_estim_indx,j_estim_indx),gamma_2(i_estim_indx,j_estim_indx),gamma_3(i_estim_indx,j_estim_indx)];
        Wei = eye(5);
%         Wei = diag([beta_1_std(i_estim_indx,j_estim_indx);beta_2_std(i_estim_indx,j_estim_indx);gamma_1_std(i_estim_indx,j_estim_indx);gamma_2_std(i_estim_indx,j_estim_indx);gamma_3_std(i_estim_indx,j_estim_indx)].^(-2));
        error(i_estim_indx,j_estim_indx) = (X_temp - X_target)*Wei*(X_temp' - X_target');
        
    end
end



[~,min_index] = min(error(:));
[min_index_row, min_index_col] = ind2sub(size(error),min_index);

estimated_rho = rho_range(min_index_row)
estimated_lambda = lambda_range(min_index_col)



% note: even in the optimal estimation, investment is highly sensitive to
% profit. a unit increase in z, allows firm to better hedge against ad-hoc
% default rule, so firm whil scale up crazily, even though it's inefficient
% in terms of operational scale, just for the sake of raising debt as well
% and enjoying the taxshield benefit.
