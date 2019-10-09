% Shasha Wang PS1 Q2 FNCE 937
clc
clear all;
close all;

cd 'E:\Dropbox\fall 19-20\Finance 937\PS1\question_2'

%% Q2 Solving the problem of the firm in continuous time
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
tolerance = 10^(-5);
% distance = 100;% I drop this line because I will not use while loop
step_size= 0.002;
%Delta = 1000;
distance = 100*ones(max_iter,1);

dVforward = zeros(Nk,1);
dVbackward = zeros(Nk,1);
dVSteadyState = ones(Nk,1); % in steady state k_dot=0, hence due to FOC k/b*(v'-1)=0, hence v'=1
investment = zeros(Nk,1);
%c = zeros(Nk,1);

%INITIAL GUESS
%v0 = (a.*k.^alpha).^(1-s)/(1-s)/r;
value0 = ones(Nk,1);
value = value0;
%iteration = 1;

tic;
%max_iter=10;
for n=1:max_iter
    V_old = value;
    % forward difference
    dVforward(1:Nk-1) = (V_old(2:Nk)-V_old(1:Nk-1))/d_k;
    %dVforward(I) = (a.*k_max.^alpha - delta.*k_max)^(-s); %state constraint, for stability
    % backward difference
    dVbackward(2:Nk) = (V_old(2:Nk)-V_old(1:Nk-1))/d_k;
    %dVbackward(1) = (Aprod.*kmin.^a - d.*kmin)^(-s); %state constraint, for stability
    
    k_dot_forward = k_grid/b .* (dVforward - 1);
    k_dot_backward = k_grid/b .* (dVbackward - 1);
    
    %investment_forward = 
        
    %consumption and savings with forward difference
    %cf = dVforward.^(-1/s);
    %muf = Aprod.*k.^a - d.*k - cf;
    %consumption and savings with backward difference
    %cb = dVbackward.^(-1/s);
    %mub = Aprod.*k.^a - d.*k - cb;
    %consumption and derivative of value function at steady state
    %c0 = Aprod.*k.^a - d.*k;
    %dV0 = c0.^(-s);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    indicator_forward = k_dot_forward > 0; %below steady state
    indicator_backward = k_dot_backward < 0; %above steady state
    indicator_center = (1-indicator_forward-indicator_backward); %at steady state

    dV_Upwind = dVforward.*indicator_forward + dVbackward.*indicator_backward + dVSteadyState.*indicator_center; %important to include third term
    %c = dV_Upwind.^(-1/s);
    %u = c.^(1-s)/(1-s);
    
    %CONSTRUCT MATRIX
    X = -min(k_dot_backward,0)/d_k;
    Y = -max(k_dot_forward,0)/d_k + min(k_dot_backward,0)/d_k;
    Z = max(k_dot_forward,0)/d_k;
    
    %full matrix: slower
    %     for i=2:Nk-1
    %         A(i,i-1) = x(i);
    %         A(i,i) = y(i);
    %         A(i,i+1) = z(i);
    %     end
    %     A(1,1)=y(1); A(1,2) = z(1);
    %     A(Nk,Nk)=y(Nk); A(Nk,Nk-1) = x(Nk);
   
    %sparse matrix: faster
    %A =spdiags(Y,0,Nk,Nk)+spdiags(X(2:Nk),-1,Nk,Nk)+spdiags([0;Z(1:Nk-1)],1,Nk,Nk);
    A =spdiags(Y,0,Nk,Nk)+spdiags([0;X(2:Nk)],-1,Nk,Nk)+spdiags([Z(1:Nk-1);0],1,Nk,Nk);

    %B = (r + 1/Delta)*speye(I) - A;
    investment = k_grid/b.*(dV_Upwind-1)+delta*k_grid;
    divident = a*k_grid.^alpha - investment - b/2*k_grid.*(investment./k_grid-delta).^2;
    value = (divident+A*V_old-r*V_old)*step_size+V_old;
    
    %b = u + V/Delta;
    %V = B\b; %SOLVE SYSTEM OF EQUATIONS
    %Vchange = V - value;
    %value = V;   

    distance(n) = max(abs(V_old - value));
    if distance(n)<tolerance
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
    %iteration = iteration+1;
    disp(n)
end
toc;






%% Graphs
set(gca,'FontSize',14)
plot(distance,'LineWidth',2)
grid
xlabel('Iteration')
ylabel('||V^{n+1} - V^n||')

kdot = Aprod.*k.^a - d.*k - c;
Verr = c.^(1-s)/(1-s) + dV_Upwind.*kdot - r.*V_old;

set(gca,'FontSize',14)
plot(k,Verr,'LineWidth',2)
grid
xlabel('k')
ylabel('Error in HJB Equation')
xlim([kmin kmax])

set(gca,'FontSize',12)
plot(k,V_old,'LineWidth',2)
grid
xlabel('k')
ylabel('V(k)')
xlim([kmin kmax])

set(gca,'FontSize',14)
plot(k,c,'LineWidth',2)
grid
xlabel('k')
ylabel('c(k)')
xlim([kmin kmax])

set(gca,'FontSize',14)
plot(k,kdot,k,zeros(1,I),'--','LineWidth',2)
grid
xlabel('$k$','FontSize',16,'interpreter','latex')
ylabel('$s(k)$','FontSize',16,'interpreter','latex')
xlim([kmin kmax])
%print -depsc HJB_NGM.eps




