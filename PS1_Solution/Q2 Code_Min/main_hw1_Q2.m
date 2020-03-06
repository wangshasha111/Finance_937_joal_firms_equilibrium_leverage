%% HW 1 Question 2

clc; clear all; close all;
addpath('./Tauchen', './output', './functions')

%%
% Declare parameter values
alpha = 0.3;
a     = 1;
b     = 0.5;
delta = 0.05;
r     = 0.05;

% Declare capital grid parameters
I     = 500;
kmin  = 0.01;
kmax  = 10;

% Create capital grid
k  = linspace(kmin,kmax,I)'; % Capital grid column vec
dk = (kmax-kmin) / (I-1);   % Step size

% Iteration parameters
maxit  = 100;
crit   = 10^(-6);
step   = 1000;

% Setup blank vectors
dV_f = zeros(I,1); % column vec
dV_b = zeros(I,1);
inv  = zeros(I,1);


% INITIAL GUESS 
v0 = a * k.^(alpha);
v = v0;




%% Main loop for Part a


% maxit = 10; % Check with lower values first
for n = 1:maxit
    V = v;
    
    % Compute forward difference of va lue function
    dV_f(1:I-1) = (V(2:I) - V(1:I-1))/dk;
    dV_f(I)     = 1; % state constraint
    
    % Compute backward difference of value function
    dV_b(2:I) = (V(2:I) - V(1:I-1))/dk;
    dV_b(1)   = 1; % state constraint
    
    % Investment and capital policy with forward difference
    inv_f    = ((dV_f - 1)/b + delta) .* k;
    kDrift_f = inv_f - delta*k;
    
    % Investment and capital policy with backward difference
    inv_b    = ((1/b) * (dV_b - 1) + delta) .* k;
    kDrift_b = inv_b - delta*k;
    
    % Steady state:
    inv_0 = delta * k;
    dV_0  = 1 + b*(inv_0 ./ k - delta);
    
    % Indicator vectors for f/b difference
    If = kDrift_f > 0;
    Ib = kDrift_b < 0;
    I0 = (1-If-Ib); % 1 if neither of the above is true
    
    % Note: dV_upwind makes a choice of forward or backward difference depending
    % on the sign of the drift.
    
    % Grab forward backward diff's for the relavant regions
    dV_upwind = dV_f.*If + dV_b.*Ib + dV_0.*I0;
    inv       = ((1/b) * (dV_upwind - 1) + delta) .* k;
    d         = a*k.^(alpha) - inv - phi(inv,k,b,delta);
    kDrift    = kDrift_f.*If + kDrift_b.*Ib + 0.*I0; % Let zeros populate where neither backward or forward applies 
    
    % X,Y,Z are block vectors that form A
    X = -min(kDrift_b,0)/dk;
    Y = -max(kDrift_f,0)/dk + min(kDrift_b,0)/dk;
    Z = max(kDrift_f,0)/dk;
    
    % Create sparse matrix A
    A = spdiags(Y,0,I,I) + spdiags(X(2:I), -1, I,I) + spdiags([0;Z(1:I-1)],1,I,I);
    
    % Setup the iterative system to be solved at each nth iteration:
    %   B(n) v(n+1) = C(n)
    % where
    %   B(n)   = 1/step + r - A(n)
    %   v(n+1) = next step value function
    %   C(n)   = d(n) + v(n) / step
    B = (1/step + r)*speye(I) - A;
    C = d + V/step;
    
    % Solve for v(n+1)
    V = B\C;
    
    % Compute update error
    err = max(abs(V-v));
    v   = V;
    
    % Show iteration
    disp(strcat('Iteration = ',num2str(n)))
    
    % Convergence criteria
    if err < crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end    
end

    
%% Plot resulting graphs for Part (a)


fig_Q2_value_func(k,V,  'output/fig_Q2_value_func.pdf')
fig_Q2_inv_policy(k,inv,'output/fig_Q2_inv_policy.pdf')
fig_Q2_kDrift(k,kDrift, 'output/fig_Q2_kDrift.pdf')


%% Part (b): Add Poisson Uncertainty

% Structure for Poisson process of a
a = [0.9; 1.0; 1.1];
P = [0.5 0.5 0; 0.25 0.5 0.25; 0 0.5 0.5];
S = 3; % number of states

% Otherwise same parameters

% Declare capital grid parameters
I     = 500;
kmin  = 0.01;
kmax  = 10;

% Create capital grid
k  = linspace(kmin,kmax,I)'; 
dk = (kmax-kmin) / (I-1);   % Step size

% Construct k grid MATRIX
k_mat = k*ones(1,3); % Each column is k(a)

% Iteration parameters
maxit  = 100;
crit   = 10^(-6);
step   = 1000;

% Setup blank vectors
dV_f = zeros(I,S); % NOW A MATRIX
dV_b = zeros(I,S);
inv  = zeros(I,S);

% Setup constant sparse matrix Aswitch. Used to get probabilities in A
% (IS)x(IS) matrix of blocks.
Aswitch = [zeros(I)        P(1,2)*speye(I) P(1,3)*speye(I); ...
           P(2,1)*speye(I) zeros(I)        P(2,3)*speye(I); ...
           P(3,1)*speye(I) P(3,2)*speye(I) zeros(I)];

% INITIAL GUESS 
v0 = k.^(alpha) * a'; % column x row .= matrix
v = v0;

%% Main loop Part (b)

% maxit = 10; % Check with lower values first
for n = 1:maxit
    V = v;
    
    % Compute forward difference of va lue function
    dV_f(1:(I-1),:) = (V(2:I,:) - V(1:(I-1),:))/dk;
    dV_f(I,:)     = 1; % state constraint for the last row
    
    % Compute backward difference of value function
    dV_b(2:I,:) = (V(2:I,:) - V(1:I-1,:))/dk;
    dV_b(1,:)   = 1; % state constraint for the first row
        
    % Investment and capital policy with forward difference
    inv_f    = ((dV_f - 1)/b + delta) .* k_mat;
    kDrift_f = inv_f - delta*k_mat;
    
    % Investment and capital policy with backward difference
    inv_b    = ((1/b) * (dV_b - 1) + delta) .* k_mat;
    kDrift_b = inv_b - delta*k_mat;
    
    % Steady state:
    inv_0 = delta * k_mat;
    dV_0  = 1 + b*(inv_0 ./ k_mat - delta);
    
    % Indicator vectors for f/b difference
    If = kDrift_f > 0;
    Ib = kDrift_b < 0;
    I0 = (1-If-Ib); % 1 if neither of the above is true
    
    % Note: dV_upwind makes a choice of forward or backward difference depending
    % on the sign of the drift.
    
    % Grab forward backward diff's for the relavant regions
    dV_upwind = dV_f.*If + dV_b.*Ib + dV_0.*I0;
    inv       = ((1/b) * (dV_upwind - 1) + delta) .* k_mat;
    d         = k.^(alpha) * a' - inv - phi(inv,k_mat,b,delta);
    kDrift    = kDrift_f.*If + kDrift_b.*Ib + 0.*I0; % Let zeros populate where neither backward or forward applies 
    
    % Sparse matrix A now looks different. Includes transition matrix P
    % Following notes, construct X,Y,Z matrices of size (IxS)
    X        = -min(kDrift_b,0)/dk;
    Z        = max(kDrift_f,0)/dk;
    
    % Constructing Y takes care. Manipulate P
    non_diag = P-diag(diag(P));            % same matrix as P with zeros on the diagonal
    P_other  = ones(I,1)*sum(non_diag,2)'; % sum of probability of exiting current state. same value across rows.
    Y        = -max(kDrift_f,0)/dk + min(kDrift_b,0)/dk - P_other;
     
    % Create block matrices for sparse matrix A
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A3=spdiags(Y(:,3),0,I,I)+spdiags(X(2:I,3),-1,I,I)+spdiags([0;Z(1:I-1,3)],1,I,I);
    
    % Create sparse matrix A using blocks and Aswitch
    A = [A1,          sparse(I,I), sparse(I,I); ...
         sparse(I,I), A2,          sparse(I,I); ...
         sparse(I,I), sparse(I,I), A3         ] + Aswitch;
    
    % Setup the linear system to be solved at each nth iteration:
    %   B(n) v(n+1) = C(n)
    % where
    %   B(n)   = (1/step + r) - A(n)
    %   v(n+1) = next step value function
    %   C(n)   = d(n) + v(n) / step
    B = (1/step + r)*speye(3*I) - A;
    
    % Stack columns to get a vector
    d_stack = d(:);
    V_stack = v(:);
    
    % Define C
    C = d_stack + V_stack/step;
    
    % Solve for v(n+1)
    V_stack = B\C;
    
    % Convert back to matrix
    V = reshape(V_stack, [I,3]); 
    
    % Compute update error
    err = max(abs(V-v));
    v   = V;
    
    % Show iteration
    disp(strcat('Iteration = ',num2str(n)))
    
    % Convergence criteria
    if err < crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end    
end

%% Plot functions

fig_Q2b_value_func(k,V,  'output/fig_Q2b_value_func.pdf')
fig_Q2b_inv_policy(k,inv,'output/fig_Q2b_inv_policy.pdf')
fig_Q2b_kDrift(k,kDrift, 'output/fig_Q2b_kDrift.pdf')


%% Misc functions

function adj_cost = phi(i,k,b,delta)
% Compute adjustment costs


adj_cost = b/2 * (i./k -delta).^2 .* k;
end

