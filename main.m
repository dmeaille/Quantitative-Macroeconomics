clear all 
close all

% Parameters
alpha = 0.4;          % Capital share in production
beta = 0.95;          % Discount factor
delta = 0.1;          % Depreciation rate
sigma = 2;            % Utility parameter
tol = 1e-6;           % Tolerance for convergence
max_iter = 1000;      % Maximum iterations
theta = 2;            % Labour disutility parameter
mu = 1;               % Labour elasticity parameter

% Options for fsolve
% (not used anymore)
opts = optimset('Diagnostics','off', 'Display','off', 'MaxIter',1000,'MaxFunEval',1000);      


nbk = 50;            % Number of element on the grid for capital


% Steady-state equations
XX = ((1/beta - (1-delta))/alpha)^(1/(alpha-1)); 
% Labour steady-state
n_ss = ((1-alpha)/theta * XX^alpha / (XX^alpha - delta*XX)^sigma)^(1/(mu+sigma)); 

% Capital stock steady-state
k_ss = XX*n_ss;                         

% Consumption stock steady-state
c_ss = k_ss^alpha * n_ss^(1-alpha) - delta*k_ss; 

% Discretize capital grid
k_min = 0.9 * k_ss;
k_max = 1.1 * k_ss;
kgrid = linspace(k_min, k_max, nbk)';

% Discretize labour grid 
nbn = nbk;
n_min = 0.5*n_ss;
n_max = 1.4*n_ss; 
ngrid = linspace(n_min, n_max, nbn)';

% Discretize consumption grid (for interpolation in VFI)
nbc = 2000;
c_min = 0.5*c_ss;
c_max = 1.8*c_ss; 
cgrid = linspace(c_min, c_max, nbc)';

% Discretizing the AR(1) Process for TFP 
rho = 0.95; 
mean_A = 1; 
std_eps = 0.007; 
N_states = 5; 
[A, PI] = rouwenhorst(N_states, mean_A, rho, std_eps); 
A = A';

% to plot the deterministic transition (just finding the index with A
% productivity shock = 1)
ind_no_shock = median(1:N_states);      

% Solution of the model: Reaction from a shock
n_rep = 50;             % number of steps of the transition paths
ind_shock = 1;          % index of shock: must be < N_states (here, worst shock) 
k_0 = k_ss*A(ind_shock);         % shock on initial capital



% To interpolate on matrices: 
function ans = interp_n(kgrid, c_policy_vfi, k_next_vfi, method, N_states)

    for i=1:N_states
        ans(1,i) = interp1(kgrid, c_policy_vfi(:,i), k_next_vfi(1,i), method);
    end 
end





%% VFI no stochastic 
% The way wa proceed is the following: 
% We will use interpolation for k(t+1), based on a current grid for c(t) that
% we allow to be much richer than the grid for k(t)

clear vi dr Kgrid Cgrid N_sol c_temp Kp interp_n u_temp
close all
 


% Preparation of the iteration grid 
% The grid won't change from iteration to iteration, so that we can compute
% it just once to speed up the computations of the VFI

[Kgrid, Cgrid] = meshgrid(kgrid, cgrid); 

% Kgrid is cst in columns, Cgrid along lines
% Result: iterating over initial k comes down to iterating over columns of
% Kgrid

% We then compute the labour that is derived from each c(t) and k(t) 
for i=1:N_states
    N_sol(:,:,i) = ((1-alpha)/theta*A(i).*Kgrid.^(alpha).*Cgrid.^(-sigma)).^(1/(mu+alpha)); 
    Kp(:,:,i) = A(i)*Kgrid.^(alpha).* N_sol(:,:,i).^(1-alpha) + (1-delta).* Kgrid - Cgrid;
    c_temp(:,:,i) = A(i)*Kgrid.^(alpha).* N_sol(:,:,i).^(1-alpha) + (1-delta).* Kgrid - Kp(:,:,i);
end 
Cgrid = c_temp;
% We have thus: 
% column of Kp: possible start of period k 
% row of Kp: possible c along with optimal n in N_sol 

% Finally, we compute the utility grid
if sigma ~= 1
    util = Cgrid.^(1-sigma)/(1-sigma) - theta*N_sol.^(1+mu)/(1+mu); 
else
    util = log(Cgrid) - theta*N_sol.^(1+mu)/(1+mu);
end

% To sum up, each column of util corresponds to a given k(t), and each row
% has a potential c(t) along with the implied n(t)
% Now, we just need to iterate over this grid


% Value function initialization
v = ones(nbk, N_states).*((1/3*kgrid.^alpha).^(1-sigma) - 1)/(1-sigma); 

% Criterium of convergence - initialization 
crit = 1;

% To store intermediate results: 
tv = zeros(nbk, N_states);
dr = zeros(nbk, N_states);

% Actual VFI computation
tic
for iter=1:max_iter
    for i=1:nbk
        for j=1:N_states
            % c(t) varies among its potential values in the grid
            % resulting k(t+1): I want for a given initial k, all the combinations with c
            
            kp = Kp(:,i,j);   % corresponding column for k(i), state j of productivity
            u = util(:,i,j);  % corresponding column for k(i), state j of productivity
            u_temp(:,j) = u;  % matrix of utility, state j of productivity
            
            % VFI interpolation of the k(t+1) that we have obtained
            vi(:,j) = interp1(kgrid,v(:,j),kp,'pchip', 'extrap');
        end

        for j=1:N_states
            % we are today at z, such that we use the corresponding row of
            % PI when computing the expected next period Value function
            [tv(i,j),dr(i,j)] = max(u_temp(:,j)+beta*vi*PI(j,:)');
        end

        % Sum up: 
        % Here, we have each column of u_temp representing a possible value
        % of A today 
        % And each value of vi representing a possible value of A tomorrow
        % We need to take the maximum, given today's A, with the
        % distribution of A tomorrow, given their transition probabilities
                
    end
    
    crit = max(max(abs(tv-v)));     % Update convergence criterion
    v=tv;                           % Update the value function

    if crit < tol
        disp(['VFI Converged in ', num2str(iter), ' iterations']);
        break;
    end

end
toc


c_policy = cgrid(dr);       % Use the indices of the optimal consumption policy to find the policy function
n_policy = ((1-alpha)/theta.*A.*kgrid.^alpha .* c_policy.^(-sigma)).^(1/(mu+alpha)); % labour policy function 
k_policy = A.*kgrid.^alpha.*n_policy.^(1-alpha) +(1-delta)*kgrid-c_policy; % k(t+1) policy function 
utility_vfi= (c_policy.^(1-sigma)-1)/(1-sigma);
v = utility_vfi/(1-beta);


% Plot the results
set(gca,'DataAspectRatio',[1,1,1])

figure (1)
plot(kgrid, c_policy, 'b-', 'LineWidth', 2);
hold on 
xlabel('Capital (k)');
ylabel('Policy Functions');
legend('Consumption Policy');
title('Ramsey Growth Model with VFI');
grid on;

c_policy_vfi = c_policy; 
k_policy_vfi = k_policy; 
n_policy_vfi = n_policy; 

figure (2)
plot(kgrid, k_policy, 'r--', 'LineWidth', 2);
hold on 
xlabel('Capital (k)');
ylabel('Policy Functions');
legend('Capital Policy');
title('Ramsey Growth Model with VFI');
grid on;


figure (3)
plot(kgrid, v, 'r--', 'LineWidth', 2);
hold on 
xlabel('Capital (k)');
ylabel('Value Functions');
legend('Value function');
title('Ramsey Growth Model with VFI');
grid on;


%% VFI Euler errors 



Euler_er=zeros(nbk,N_states);
for t=1:nbk
    k_current = kgrid(t, :);
    c_vfi = c_policy_vfi(t, :);
    k_next_vfi = k_policy_vfi(t, :);
    c_next_vfi = interp_n(kgrid, c_policy_vfi, k_next_vfi, 'spline', N_states);
    n_next_vfi = interp_n(kgrid, n_policy_vfi, k_next_vfi, 'spline', N_states);
    
    for j=1:N_states
        euler_c(:,j) = (beta * (1 - delta + A .* alpha .* k_next_vfi.^(alpha - 1) .* n_next_vfi.^(1-alpha) ).* c_next_vfi.^(-sigma)).^(-1/sigma)*PI(j,:)';
    end

    Euler_er(t, :) = (euler_c-c_vfi)./c_vfi .* 100;
end

figure (4);
plot(kgrid, Euler_er, '-', 'LineWidth', 2);
hold on 
%plot(kgrid, c_policy, 'b-', 'LineWidth',1.5)
xlabel('Capital (k)');
ylabel('Euler Errors');
title('Euler errors with VFI method');
grid on;






%% STOCHASTIC EGM 

% For EGM, we need a grid for today's capital, one for tomorrow (they are
% the same here). 
k_today = kgrid; 
k_tomorow = kgrid;
% We derive a grid for today's labor from today's capital
n_today = ((1-alpha)/theta*k_today).^(1/(mu+alpha)); 


% Initial policy guess for consumption 
% note that this is a function of today's capital: this is the optimal
% consumption today, given today's capital
c_policy = A .* k_today.^alpha + (1-delta)*k_today;

tic
% Value Function Iteration with Endogenous Grid Method
for iter = 1:max_iter
    
    % First, get c tomorrow, interpolating c_policy function of today's
    % capital, with tomorrow's capital 
    c_next = interp1(k_today, c_policy, k_tomorow, 'spline', 'extrap');
    
    % Implied n(t+1) labor for tomorrow, knowing c(t+1) from just before and using
    % the grid for k(t+1)
    n_tomorow = ((1-alpha)/theta*A.*k_tomorow.^alpha .* c_next.^(-sigma)).^(1/(mu+alpha)); 

    % Compute consumption today from stochastic Euler Equation, iterating on k(t+1)
    for j=1:N_states
        c_today(:,j) = c_next .* (beta * (alpha*A.*k_tomorow.^(alpha-1).*n_tomorow.^(1-alpha) + 1 - delta)).^(-1/sigma)*PI(j,:)';
    end

    cash = c_today + k_tomorow;

    % Compute capital today from the grid we get on consumption today c(t)
    % This is our endogenous grid. We have it by solving the non-linear
    % resource constraint for k(t), knowing c(t) and k(t+1)
    for i=1:nbk
        for j=1:N_states
        k_tild(i,j) = newton(k_today(i),c_today(i,j),A(j),alpha,theta,mu,delta,sigma,cash(i,j)); 
        end
    end
    
  
    % Interpolate policy onto original grid: using the endogenous grid we
    % found, we can update our policy function for consumption 
    for j=1:N_states
        c_policy_new(:, j) = interp1(k_tild(:,j), c_today(:,j), k_today , 'spline', 'extrap');
    end
    
    % Convergence check
    if max(abs(c_policy_new - c_policy)) < tol
        disp(['EGM Converged in ', num2str(iter), ' iterations']);
        break;
    end
    
    % Update consumption policy
    c_policy = c_policy_new;
    c_policy = max(c_policy, 0);
end
toc

% Labour policy function 
n_policy = ((1-alpha)/theta*A.*k_today.^alpha .* c_policy.^(-sigma)).^(1/(mu+alpha)); 

% Capital policy function
k_policy = A.* k_today.^alpha.*n_policy.^(1-alpha) + (1-delta).*k_today - c_policy;




% Plot the policy functions
figure (5);
plot(kgrid, c_policy, 'b-', 'LineWidth', 2);
xlabel('Capital (k)');
ylabel('Policy Functions');
legend('Consumption Policy');
title('Ramsey Growth Model with Endogenous Grid Method');
grid on;
% 
figure (6);
plot(kgrid, k_policy, 'r--', 'LineWidth', 2);
xlabel('Capital (k)');
ylabel('Policy Functions');
legend('Capital Policy');
title('Ramsey Growth Model with Endogenous Grid Method');
grid on;

c_policy_egm = c_policy; 
k_policy_egm = k_policy; 
n_policy_egm = n_policy;



%% Euler errors



Euler_er=zeros(nbk,N_states);
for t=1:nbk
    k_current = k_today(t);
    c_egm = c_policy_egm(t, :);
    k_next_egm = k_policy_egm(t, :);
    c_next_egm = interp_n(kgrid, c_policy_egm, k_next_egm, 'spline', N_states);
    n_next_egm = interp_n(kgrid, n_policy_egm, k_next_egm, 'spline', N_states);
    
    for j=1:N_states
        euler_c(:,j) = (beta * (1 - delta + A .* alpha .* k_next_egm.^(alpha - 1) .* n_next_egm.^(1-alpha) ).* c_next_egm.^(-sigma)).^(-1/sigma)*PI(j,:)';
    end

    Euler_er(t, :) = (euler_c-c_egm)./c_egm .* 100;
end

figure (7);
plot(kgrid, Euler_er, '-', 'LineWidth', 2);
hold on 
%plot(kgrid, c_policy, 'b-', 'LineWidth',1.5)
xlabel('Capital (k)');
ylabel('Euler Errors, in %');
title('Euler errors with EGM method');
grid on;



%% Analysis 

disp('Checking that we get similar results from the two methods for the policy function for consumption:')
disp(c_policy_vfi-c_policy_egm)



%% SOLUTION BY LOG-LINEARIZATION 


% Output steady-state
y_ss=k_ss^alpha*n_ss^(1-alpha);
% Marginal Productivity
margprod = alpha*k_ss^(alpha-1)*n_ss^(1-alpha);

a= [[beta*(alpha-1)*margprod,beta*margprod,-sigma,beta*(1-alpha)*margprod]; 
    [0, 0,0,0]; 
    [k_ss,0,0,0];
    [0,1,0,0]];

b= [[0,0,-sigma,0]; 
    [alpha, 1, -sigma, -(mu+alpha)]; 
    [k_ss/beta,y_ss,-c_ss,(1-alpha)*y_ss];
    [0,rho,0,0]];


nk=2; %number of state variable
[f,p] = solab(a,b,nk);
% We extract consumption and labor policies
c_polfunc=f(1,:);
n_polfunc=f(2,:);
% We extract the law of motion for capital
LOM=p(1,:); 


% To make it easy to compare with VFI and EGM, we compute the equivalent
% policy functions: 
c_policy_loglin = c_ss*exp(c_polfunc(1)*log(kgrid/k_ss)+c_polfunc(2)*log(A/1));
n_policy_loglin = n_ss*exp(n_polfunc(1)*log(kgrid/k_ss)+n_polfunc(2)*log(A/1));
k_policy_loglin = k_ss*exp(LOM(1)*log(kgrid/k_ss)+LOM(2)*log(A/1));





% Plot the policy functions
figure (8);
plot(kgrid, c_policy_loglin, 'b-', 'LineWidth', 2);
xlabel('Capital (k)');
ylabel('Policy Functions');
legend('Consumption Policy');
title('Ramsey Growth Model with Log-linearization');
grid on;
% 
figure (9);
plot(kgrid, k_policy_loglin, 'r--', 'LineWidth', 2);
xlabel('Capital (k)');
ylabel('Policy Functions');
legend('Capital Policy');
title('Ramsey Growth Model with Log-linearization');
grid on;




%% PLOTS OF THE REACTIONS TO A SHOCK 


k_transit = zeros(n_rep, 3);        % initialize the path vector of capital
c_transit = zeros(n_rep, 3);        % initialize the path vector of consumption
n_transit = zeros(n_rep, 3);        % initialize the path vector of labor



% Take the closest value to k_0 we have in our grid, interpolating to nearest
k_init = interp1(kgrid, kgrid, k_0, 'nearest'); 


for i=1:n_rep
    if i==1         
        k_transit(1, 1)=k_init; % initialize the transition path of capital
        temp=k_transit(i, 1);
    else
        % start from the initial value and apply the policy function
        % that is, interpolate with the previous value for capital
        k_transit(i, 1)=interp1(kgrid, k_policy_loglin(:,ind_no_shock), temp, 'spline');
        temp=k_transit(i, 1);      
    end
    % use policy function for consumption and labor
    n_transit(i,1) = interp1(kgrid, n_policy_loglin(:,ind_no_shock), temp, 'spline'); 
    c_transit(i,1) = interp1(kgrid, c_policy_loglin(:,ind_no_shock), temp, 'spline');
end



% VFI                     
for i=1:n_rep
    if i==1         
        k_transit(1, 2)=k_init;
        temp=k_transit(i, 2);
    else
        k_transit(i, 2)=interp1(kgrid, k_policy_vfi(:,ind_no_shock), temp, 'spline');
        temp=k_transit(i, 2);      
    end

    n_transit(i,2) = interp1(kgrid, n_policy_vfi(:,ind_no_shock), temp, 'spline'); 
    c_transit(i,2) = interp1(kgrid, c_policy_vfi(:,ind_no_shock), temp, 'spline');
end


% EGM
for i=1:n_rep
    if i==1         
        k_transit(1, 3)=k_init;
        temp=k_transit(i, 3);
    else
        k_transit(i, 3)=interp1(kgrid, k_policy_egm(:,ind_no_shock), temp, 'spline', 'extrap');
        temp=k_transit(i, 3);      
    end

    n_transit(i,3) = interp1(kgrid, n_policy_egm(:,ind_no_shock), temp, 'spline', 'extrap'); 
    c_transit(i,3) = interp1(kgrid, c_policy_egm(:,ind_no_shock), temp, 'spline', 'extrap');
end



%%% Plotting the resulting transition paths

figure (10)
% VFI 
subplot(1, 3, 1)
yyaxis left
plot(1:n_rep, c_transit(1:n_rep,1))
hold on 
yline(c_ss, '--')
yyaxis right
plot(1:n_rep, k_transit(1:n_rep, 1))
yline(k_ss, '--')
title('Log-lin Model')
hold off
% EGM
subplot(1, 3, 2)
yyaxis left
plot(1:n_rep, c_transit(1:n_rep,2))
hold on 
yline(c_ss, '--')
yyaxis right
plot(1:n_rep, k_transit(1:n_rep,2))
yline(k_ss, '--')
xlabel('time')
title('VFI Model')
hold off
% Log-lin
subplot(1, 3, 3)
yyaxis left
plot(1:n_rep, c_transit(1:n_rep,3))
hold on 
yline(c_ss, '--')
yyaxis right
plot(1:n_rep, k_transit(1:n_rep,3))
yline(k_ss, '--')
title('EGM Model')
legend('Consumption', 'Capital', 'Location', 'northeast')
hold off






%% SIMULATION OF MARKOV-CHAINS 



T = 1000;           % Number of periods 
N_sim = 100;        % Number of simulations
rng(0) % fix state of random numer generator

% Continuous Markov Chain for reference, we won't use it
a_cont = zeros(T, N_sim); 
for j=1:N_sim
    a_cont(j,1) = 0;
    for t=2:T
        a_cont(j,t) = rho*a_cont(j,t-1) + normrnd(0,std_eps);
    end
end


% Discrete Markov Chain
a_disc = zeros(T, N_sim);   % one simulation is one column 
a_disc(1,:) = ind_no_shock; % We start with A = 1, that is no shock (or steady state) for all chains

for t=2:T
    for sim=1:N_sim
        a_disc(t,sim) = randsample(1:N_states, 1, true, PI(a_disc(t-1, sim),:)); % We weight by the probabilities coming from last period state
    end
end

% Equivalent path for productivity: 
A_disc = A(a_disc);


%% PLOT OF THE SIMULATIONS 

k_sim_vfi = zeros(T, N_sim);
k_sim_egm = zeros(T, N_sim);

c_sim_vfi = zeros(T, N_sim);
c_sim_egm = zeros(T, N_sim);

n_sim_vfi = zeros(T, N_sim);
n_sim_egm = zeros(T, N_sim);


k_sim_vfi(1,:) = k_ss;
k_sim_egm(1,:) = k_ss;

for t=1:T
    for sim=1:N_sim
        shock_ind = a_disc(t,sim);

        % We simply interpolate, based on the index of the shock: 
        % Since we are in t in a state for A, the probabilities for next
        % period distribution will vary accordingly: P(A'|A) 
        % Hence, we need to choose the correct policy function, since they
        % depend on today's state
        
        if t < T
            k_sim_vfi(t+1, sim) = interp1(kgrid, k_policy_vfi(:,shock_ind), k_sim_vfi(t, sim), 'linear', 'extrap');
            k_sim_egm(t+1, sim) = interp1(kgrid, k_policy_egm(:,shock_ind), k_sim_egm(t, sim), 'linear', 'extrap');
        end
        
        c_sim_vfi(t, sim) = interp1(kgrid, c_policy_vfi(:,shock_ind), k_sim_vfi(t, sim), 'linear', 'extrap');
        c_sim_egm(t, sim) = interp1(kgrid, c_policy_egm(:,shock_ind), k_sim_egm(t, sim), 'linear', 'extrap');

        n_sim_vfi(t, sim) = interp1(kgrid, n_policy_vfi(:,shock_ind), k_sim_vfi(t, sim), 'linear', 'extrap');
        n_sim_egm(t, sim) = interp1(kgrid, n_policy_egm(:,shock_ind), k_sim_egm(t, sim), 'linear', 'extrap');

    end
end


% Choose a simulation randomly
nb_sim = 5; 
figure (11)

% K 
subplot(2, 1, 1)
yyaxis left
plot(1:T, k_sim_vfi(:,nb_sim), 'r-'); 
hold on 
plot(1:T, k_sim_egm(:,nb_sim), 'b-'); 
yline(k_ss, '--')
yyaxis right
plot(1:T, A_disc(:, nb_sim));
title('Simulated Capital path');
legend('VFI', 'EGM', 'Shock');
hold off
% C
subplot(2, 1, 2)
plot(1:T, c_sim_vfi(:,nb_sim)); 
hold on 
plot(1:T, c_sim_egm(:,nb_sim)); 
yline(c_ss, '--')
title('Simulated Consumption path')
hold off
% N
% subplot(1, 3, 3)
% plot(1:T, n_sim_vfi(:,nb_sim)); 
% hold on 
% plot(1:T, n_sim_egm(:,nb_sim)); 
% yline(n_ss, '--')
% title('Simulated Labour path')
% hold off




%% ESTIMATION OF THE MOMENTS

% Standard deviation of output
output_sim_vfi = A_disc.*k_sim_vfi.^(alpha) .* n_sim_vfi.^(1-alpha); 
output_sim_egm = A_disc.*k_sim_egm.^(alpha) .* n_sim_egm.^(1-alpha);

std_output_vfi = mean(std(output_sim_vfi(100:T,:)));
std_output_egm = mean(std(output_sim_egm(100:T,:)));

% Standard deviation of consumption
std_c_vfi = mean(std(c_sim_vfi(100:T,:)));
std_c_egm = mean(std(c_sim_egm(100:T,:)));

% Standard deviation of capital
std_k_vfi = mean(std(k_sim_vfi(100:T,:)));
std_k_egm = mean(std(k_sim_egm(100:T,:)));

% Standard deviation of labor
std_n_vfi = mean(std(n_sim_vfi(100:T,:)));
std_n_egm = mean(std(n_sim_egm(100:T,:)));



% Correlation with output 
% consumption
for sim=1:N_sim
    correl_output_c_vfi(sim) = corr(output_sim_vfi(:,sim), c_sim_vfi(:,sim));
    correl_output_c_egm(sim) = corr(output_sim_egm(:,sim), c_sim_egm(:,sim));

end
correl_output_consumption_vfi = mean(correl_output_c_vfi); 
correl_output_consumption_egm = mean(correl_output_c_egm); 


% capital
for sim=1:N_sim
    correl_output_k_vfi(sim) = corr(output_sim_vfi(:,sim), k_sim_vfi(:,sim));
    correl_output_k_egm(sim) = corr(output_sim_egm(:,sim), k_sim_egm(:,sim));

end
correl_output_capital_vfi = mean(correl_output_k_vfi); 
correl_output_capital_egm = mean(correl_output_k_egm); 


% labor
for sim=1:N_sim
    correl_output_n_vfi(sim) = corr(output_sim_vfi(:,sim), n_sim_vfi(:,sim));
    correl_output_n_egm(sim) = corr(output_sim_egm(:,sim), n_sim_egm(:,sim));

end
correl_output_labor_vfi = mean(correl_output_n_vfi); 
correl_output_labor_egm = mean(correl_output_n_egm); 



%% Euler error on the simulated path


Euler_er=zeros(N_sim-1,2);
for t=2:N_sim
    %left hand side VFI and EGM:
    lhs_vfi=c_sim_vfi(nb_sim,t-1)^(-sigma);
    lhs_egm=c_sim_egm(nb_sim,t-1)^(-sigma);
    %right hand side VFI and EGM:
    rhs_vfi = beta * (1 - delta + A_disc(nb_sim,t) * alpha * k_sim_vfi(nb_sim,t)^(alpha - 1) .* n_sim_vfi(nb_sim,t)^(1 - alpha))* c_sim_vfi(nb_sim,t)^(-sigma);
    rhs_egm = beta * (1 - delta + A_disc(nb_sim,t) * alpha * k_sim_egm(nb_sim,t)^(alpha - 1) .* n_sim_egm(nb_sim,t)^(1 - alpha))* c_sim_egm(nb_sim,t)^(-sigma);
    Euler_er(t,1)=lhs_vfi-rhs_vfi; %error in VFI
    Euler_er(t,2)=lhs_egm-rhs_egm; %error in EGM
end

figure (12)
plot(Euler_er(:,1), 'b-'); 
hold on 
plot(Euler_er(:,2), 'r-'); 
title("Euler Error on the simulated path of consumption for VFI and EGM")
legend("VFI", "EGM")
hold off




