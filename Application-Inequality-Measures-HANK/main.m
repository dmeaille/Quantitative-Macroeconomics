clc
clear all
close all

% Setting paths
addpath '/MATLAB Drive/Quant Macro/Macro simu'
cd '/MATLAB Drive/Quant Macro/Macro simu'


% Loading the parameters
par = parameters(); 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               Household Bellman Problem using EGM              %%%                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we initialize the value function
V = (par.Amesh + par.y').^par.eis;

% We choose a starting r
r = (1/par.beta-1)*0.95; % to ensure no infinite holding of assets


% Using EGM to solve the problem 
[A, ~, a_policy, c_policy, ~] = household_ss(r, V, par);


figure
subplot(1, 2, 1)
plot(par.agrid, a_policy, '-', 'LineWidth', 0.8); 
xlabel('Asset');
ylabel('Savings');
title('Savings Policy')
% legend('poor', '', 'middle', '', 'high');
grid on
subplot(1, 2, 2)
plot(par.agrid, c_policy, '-', 'LineWidth', 0.8);
xlabel('Asset');
ylabel('Consumption');
legend('$\epsilon_g$', '$\epsilon_b$', '$\epsilon_u$','Interpreter','latex');
title('Consumption Policy');
grid on;
sgtitle("Policy Functions Using EGM")



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 General Equilibrium for r and w                %%%                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Firm labor demand depends on wage level
w = par.Z;         % wage level
N = par.labor;     % firm aggregate labour demand
Gamma = par.Gamma; % spread of mu_g on mu_b

% Error metrics
diff_metric = 1; 

% Initializing bounds for interest rate
r_do = 0;
r_up = 0.1;

% Initializing value function
V = (par.Amesh + par.y').^par.eis;


% Bisection algorithm on r
while diff_metric > 2e-7
    % new test point
    r = (r_do + r_up)/2; 

    % induced tax rate: governement budget constraint
    K = ((1+r)/par.alpha)^(1/(1-par.alpha))*N; 
    w = (1-par.alpha) * par.Z * K^(par.alpha) * N^(-par.alpha);
    update_grids 

    % EGM and computation of the distribution 
    [A, dist, a_policy, c_policy, V] = household_ss(r, V, par);

    % updating the bounds
    if A > K
        r_up = r;
    else
        r_do = r;
    end
    
    % updating the difference
    diff_metric = abs(A-K);
    % disp(diff_metric);

end

par.r = r; 
par.wage = w; 

data = table([A; K; N; r; w], 'VariableNames', {'Equilibrium Values'}, 'RowNames', {'A', 'K', 'N', 'r', 'w'});

% Display the table
disp(data);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               Stationary distribution of wealth                %%%                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plotting the Asset Distribution 
distribution_by_states = reshape(dist, par.nb_states, par.nba);
distribution_global = sum(distribution_by_states); 



figure
plot(par.agrid, distribution_global, LineWidth=1, LineStyle="-", Color="blue")
xlim([par.Amin, par.Amax])
title("Steady-state Distribution of Asset holdings")
xlabel("Assets")
ylabel("Share of Total Population")




figure 
plot(par.agrid, distribution_by_states(1,:), LineWidth=1, LineStyle="-")
hold on 
plot(par.agrid, distribution_by_states(2,:), LineWidth=1, LineStyle="-")
plot(par.agrid, distribution_by_states(3,:), LineWidth=1, LineStyle="-")
title("Distribution by Employment Status")
legend('Employed', 'Bad', 'Unemployed');



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        Measuring Inequality                    %%%                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the resulting table 
clear data
data = table([0;0;0], 'VariableNames', {'Inequality Index'}, 'RowNames', {'Gini', 'Asset top10/btm10', 'Consumption top10/btm10'});



Gamma = par.Gamma; 

Gam = linspace(0, 1, 6);

for G = Gam 
    Gamma = G; 
    update_grids

    [gini, asset_ratio, consumption_ratio] = compute_inequality_metrics(par);
    new_column = [gini; asset_ratio; consumption_ratio];

    %%% Updating the table 
    G_formatted = sprintf('%.1f', G); % Format G to two decimal places
    data = [data table(new_column, 'VariableNames', {sprintf('Gamma = %s', G_formatted)})];

end 

data(:, 1) = [];

% Display the table
disp(data);





