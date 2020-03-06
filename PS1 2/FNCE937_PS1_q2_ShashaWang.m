% Shasha Wang PS1 Q2 FNCE 937
clc
clear all;
close all;

cd 'E:\Dropbox\fall 19-20\Finance 937\PS1\question_2'

%% Q2 Solving the problem of the firm in continuous time - implicit method
% Written by Benjamin Moll
% Based on code by Benjamin Moll, using theory from
% Convergent Difference Schemes for Nonlinear Elliptic and Parabolic Equations

%% Parametization
alpha = 0.3;
b = 0.5;
delta = 0.05;
r = 0.05;
a = 1;

% k_grid
Nk = 200;
% center should be set at
k_SteadyState = (alpha*a/(r+delta))^(1/(1-alpha)); % i.e., steady state: F'(k)=r+delta
k_max = (a/delta)^(1/(1-alpha));
k_max = min(2*k_SteadyState, k_max); % Tighten the grid
k_min = 0.001*k_SteadyState;
k_grid = linspace(k_min,k_max,Nk)';
d_k = (k_max-k_min)/(Nk-1);

max_iter=10000;
tolerance = 10^(-6);
% distance = 100;% I drop this line because I will not use while loop
step_size= 0.02;
%Delta = 1000;
distance = zeros(max_iter,1);

dVforward = zeros(Nk,1);
dVbackward = zeros(Nk,1);
dVSteadyState = ones(Nk,1); % in steady state k_dot=0, hence due to FOC k/b*(v'-1)=0, hence v'=1
%investment = zeros(Nk,1);
%c = zeros(Nk,1);

%INITIAL GUESS
value0 = [1:Nk].^0.5';
value = value0;

tic;
for n=1:max_iter
    V_old = value;
    % forward difference
    dVforward(1:Nk-1) = (V_old(2:Nk)-V_old(1:Nk-1))/d_k;
    % backward difference
    dVbackward(2:Nk) = (V_old(2:Nk)-V_old(1:Nk-1))/d_k;
    
    k_dot_forward = k_grid/b .* (dVforward - 1);
    k_dot_backward = k_grid/b .* (dVbackward - 1);
    
    indicator_forward = k_dot_forward > 0; %below steady state
    indicator_backward = k_dot_backward < 0; %above steady state
    indicator_center = (1 - indicator_forward - indicator_backward); %at steady state

    dV_Upwind = dVforward.*indicator_forward + dVbackward.*indicator_backward + dVSteadyState.*indicator_center; 
    
    %CONSTRUCT MATRIX
    X = -min(k_dot_backward,0)/d_k;
    Y = -max(k_dot_forward,0)/d_k + min(k_dot_backward,0)/d_k;
    Z = max(k_dot_forward,0)/d_k;
    
    %AA = zeros(Nk,Nk);
    
    %full matrix: slower
%          for i=2:Nk-1
%              AA(i,i-1) = X(i);
%              AA(i,i) = Y(i);
%              AA(i,i+1) = Z(i);
%          end
%          AA(1,1)=Y(1); AA(1,2) = Z(1);
%          AA(Nk,Nk)=Y(Nk); AA(Nk,Nk-1) = X(Nk);
   
    %sparse matrix: faster
    A =spdiags(Y,0,Nk,Nk)+spdiags([X(2:Nk)],-1,Nk,Nk)+spdiags([0;Z(1:Nk-1)],1,Nk,Nk);

    B = (r + 1/step_size)*speye(Nk) - A;
    
    investment = k_grid/b.*(dV_Upwind-1)+delta*k_grid;
    divident = a*k_grid.^alpha - investment - b/2*k_grid.*(investment./k_grid-delta).^2;
    %value = (divident + A*V_old - r*V_old) * step_size + V_old; %explicity method
    bb = divident + V_old/step_size;
    V_new = B\bb; %SOLVE SYSTEM OF EQUATIONS
    Vchange = V_new - V_old;
    value = V_new;
    
    distance(n) = max(abs(Vchange));
    if distance(n)<tolerance
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
    disp(n)
end
toc;

%% Graphs
figure(1);
subplot(2,1,1)
set(gca,'FontSize',14)
plot(distance,'LineWidth',2)
grid
xlabel('Iteration')
xlim([0,n])
ylabel('||V^{n+1} - V^n||')
title('Value Function Convergence')
%savefig('q2a_distance')

subplot(2,1,2)
set(gca,'FontSize',14)
plot(k_grid,V_new,'LineWidth',2)
xlabel('Current Capital Stock $k$','interpreter','latex')
xlim([0,k_max])
ylabel('Value Function')
title('Value')

savefig('q2a_value_function_iteration')

figure(2)
kdot = investment - delta*k_grid;
set(gca,'FontSize',14)
plot(k_grid, kdot,k_grid, investment, 'LineWidth',2)
xlabel('Current Capital Stock $k$','interpreter','latex')
xlim([0,k_max])
ylabel('$\dot{k}$ $\:$ $i$','interpreter','latex')
title('Capital Change - $\dot{k}$ Investment - $i$','interpreter','latex')
yline(0);
xline(k_SteadyState); % Add line that makes it obvious that k_dot at k_SteadyState is 0
legend('$\dot{k}$','$i$','interpreter','latex')
savefig('q2a_kdot_investment')

%% Poisson uncertainty

% a_grid
Na=3;
a_grid=[0.9, 1, 1.1]';
a_prob=[0.5,0.5,0;0.25,0.5,0.25;0,0.5,0.5];

% k_grid
Nk = 30;
% center should be set at
k_SteadyState = (alpha*a/(r+delta))^(1/(1-alpha)); % i.e., steady state: F'(k)=r+delta
k_max = (a/delta)^(1/(1-alpha));
k_max = min(2*k_SteadyState, k_max); % Tighten the grid
k_min = 0.001*k_SteadyState;
k_grid = linspace(k_min,k_max,Nk)';
d_k = (k_max-k_min)/(Nk-1);

max_iter=20000;
tolerance = 10^(-6);
step_size= 0.02;
%Delta = 1000;
% distance = 100;% I drop this line because I will not use while loop
distance = zeros(max_iter,1);

% create transition matrix
P = kron(a_prob,speye(Nk)); % kronecker product of probability and identity matrix

% Initial Guess for value function
value0 = ones(Nk*Na,1);
value = value0;

dVforward = zeros(Nk*Na,1);
dVbackward = zeros(Nk*Na,1);
dVSteadyState = ones(Nk*Na,1); % in steady state k_dot=0, hence due to FOC k/b*(v'-1)=0, hence v'=1

investment = zeros(Nk*Na,1);
divident = zeros(Nk*Na,1);

% create sparse matrix A
A_3d = zeros(Nk,Nk,Na); % Create a 3d matrix and then stack it in a diagonal manner
A_2d = zeros(Nk*Na,Nk*Na);

tic;
for n=1:max_iter
    V_old = value;
    for ai=1:Na
        % forward difference
        dVforward((ai-1)*Nk+1:Nk*ai-1) = (V_old((ai-1)*Nk+2:Nk*ai)-V_old((ai-1)*Nk+1:Nk*ai-1))/d_k;
        % backward difference
        dVbackward((ai-1)*Nk+2:Nk*ai) = (V_old((ai-1)*Nk+2:Nk*ai)-V_old((ai-1)*Nk+1:Nk*ai-1))/d_k;

        k_dot_forward = k_grid/b .* (dVforward((ai-1)*Nk+1:Nk*ai) - 1);
        k_dot_backward = k_grid/b .* (dVbackward((ai-1)*Nk+1:Nk*ai) - 1);

        indicator_forward = k_dot_forward > 0; %below steady state
        indicator_backward = k_dot_backward < 0; %above steady state
        indicator_center = (1 - indicator_forward - indicator_backward); %at steady state

        dV_Upwind = dVforward((ai-1)*Nk+1:Nk*ai).*indicator_forward + dVbackward((ai-1)*Nk+1:Nk*ai).*indicator_backward + dVSteadyState((ai-1)*Nk+1:Nk*ai).*indicator_center; 

        %CONSTRUCT MATRIX
        X = -min(k_dot_backward,0)/d_k;
        Y = -max(k_dot_forward,0)/d_k + min(k_dot_backward,0)/d_k;
        Z = max(k_dot_forward,0)/d_k;

        %AA = zeros(Nk,Nk);

        %full matrix: slower
    %          for i=2:Nk-1
    %              AA(i,i-1) = X(i);
    %              AA(i,i) = Y(i);
    %              AA(i,i+1) = Z(i);
    %          end
    %          AA(1,1)=Y(1); AA(1,2) = Z(1);
    %          AA(Nk,Nk)=Y(Nk); AA(Nk,Nk-1) = X(Nk);

        %sparse matrix: faster
        A = spdiags(Y,0,Nk,Nk)+spdiags([X(2:Nk)],-1,Nk,Nk)+spdiags([0;Z(1:Nk-1)],1,Nk,Nk);
        A_3d(:,:,ai) = A;
        A_2d((ai-1)*Nk+1:ai*Nk,(ai-1)*Nk+1:ai*Nk) = A;
        
        investment((ai-1)*Nk+1:ai*Nk,1) = k_grid/b.*(dV_Upwind-1)+delta*k_grid;
        divident((ai-1)*Nk+1:ai*Nk,1) = a_grid(ai)*k_grid.^alpha - investment((ai-1)*Nk+1:ai*Nk,1) - b/2*k_grid.*(investment((ai-1)*Nk+1:ai*Nk,1)./k_grid-delta).^2;
        
    end
    
    % A manual way of putting numbers in the big A matrix of 2 dimensions
    % A_2d
    %A1=A_3d(:,:,1);
    %A2=A_3d(:,:,2);
    %A3=A_3d(:,:,3);
    %A_2d=[[A1,zeros(Nk,Nk),zeros(Nk,Nk)];[zeros(Nk,Nk),A2,zeros(Nk,Nk)];[zeros(Nk,Nk),zeros(Nk,Nk),A3]];
    
    A_poisson = A_2d + P - speye(Nk*Na);
    B = (r + 1/step_size)*speye(Nk*Na) - A_poisson;
   
    %investment = k_grid/b.*(dV_Upwind-1)+delta*k_grid;
    %divident = a*k_grid.^alpha - investment - b/2*k_grid.*(investment./k_grid-delta).^2;
    %value = (divident + A*V_old - r*V_old) * step_size + V_old; %explicity method
    bb = divident + V_old/step_size;
    V_new = B\bb; %SOLVE SYSTEM OF EQUATIONS
    Vchange = V_new - V_old;
    value = V_new;
    
    distance(n) = max(abs(Vchange));
    if distance(n)<tolerance
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
    disp(n)

end
toc;

%% Graphs
figure(3);
subplot(2,1,1)
set(gca,'FontSize',14)
plot(distance,'LineWidth',2)
grid
xlabel('Iteration')
xlim([0,n])
ylabel('||V^{n+1} - V^n||')
title('Value Function Convergence')
%savefig('q2b_distance')

subplot(2,1,2)
value_Nk_By_Na=reshape(value,[Nk,Na]);% turn vector to a Nk by Na matrix
set(gca,'FontSize',14)
plot(k_grid,value_Nk_By_Na,'LineWidth',0.8)
xlabel('Current Capital Stock $k$','interpreter','latex')
xlim([0,k_max])
ylabel('Value Function')
legend("low productivity","medium productivity","high productivity",'Location','southeast')
title('Value')

savefig('q2b_value_function_iteration')

figure(4)
kdot = investment - delta*repmat(k_grid,Na,1); % expand the k_grid by replication
kdot_Nk_By_Na=reshape(kdot,[Nk,Na]);
investment_Nk_By_Na=reshape(investment,[Nk,Na]);

subplot(2,1,1)
set(gca,'FontSize',14)
plot(k_grid, kdot_Nk_By_Na, 'LineWidth',0.8)
xlabel('Current Capital Stock $k$','interpreter','latex')
xlim([0,k_max])
ylabel('$\dot{k}$ $\:$ $i$','interpreter','latex')
title('Capital Change - $\dot{k}$','interpreter','latex')
yline(0);
xline(k_SteadyState); % Add line that makes it obvious that k_dot at k_SteadyState is 0
legend("low productivity","medium productivity","high productivity",'Location','southwest')
%savefig('q2a_kdot_investment')

subplot(2,1,2)
set(gca,'FontSize',14)
plot(k_grid, investment_Nk_By_Na, 'LineWidth',0.8)
xlabel('Current Capital Stock $k$','interpreter','latex')
xlim([0,k_max])
ylabel('$i$','interpreter','latex')
title('Investment - $i$','interpreter','latex')
%yline(0);
%xline(k_SteadyState); % Add line that makes it obvious that k_dot at k_SteadyState is 0
legend("low productivity","medium productivity","high productivity",'Location','southwest')
savefig('q2b_kdot_investment')