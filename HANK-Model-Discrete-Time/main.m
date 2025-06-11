clear all

% Loading the parameters
par = parameters(); 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     Steady-state distribution                  %%%                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Household aggregate labor supply 
L_hh = 1; 

% Governement debt (also asset demand)
B = par.debt; 

% Firm labor demand depends on wage level
w = par.Z;         % wage 
N = par.labor;     % firm labour demand


% Error metrics
diff_metric = 1; 

% Initializing bounds for interest rate
r_do = 0;
r_up = 0.02;

% Initializing value function
V = (par.AA + par.s').^par.eis;


% Bisection algorithm on r
tic
while diff_metric > par.tol_egm
    % new test point
    r = (r_do + r_up)/2; 

    % induced tax rate: governement budget constraint
    tax = (B*r)/(par.Z*N);

    % EGM and computation of the distribution 
    [A, dist, Anext, C, V] = household_ss(r, w, tax, V, N, par);

    % updating the bounds
    if A > B
        r_up = r;
    else
        r_do = r;
    end
    
    % updating the difference
    diff_metric = abs(A-B);
    % disp(diff_metric);

end
toc
disp("EGM problem converged to A=B")

par.r = r; 
par.tax = tax;

% C aggregate from distribution
C_agg =dist'*C(:);

% C aggregate from Aggregate resource constraint
C_hh = (1-tax) * w * L_hh + (1+r)*B - A; 
phi = 1/par.mu * (1-tax) * w * C_hh^(-par.gamma); 

% Plotting the Asset Distribution 
distribution = reshape(dist, par.M, par.N);
distribution = sum(distribution); 

plot(par.Agrid, distribution, LineWidth=1, LineStyle="-", Color="blue")
xlim([0, 20])
title("Steady-state Distribution of Asset holdings")
xlabel("Assets")
ylabel("Share of Total Population")


[Anext, C, V] = egm(r, w, tax, V, N, par);

par.phi = 1/par.mu * (1-tax) * par.Z * C_agg^(-par.gamma);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 One time Monetary Policy Shock                 %%%                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%% Computation for the steady-state reference --> should find A^hh(t) = 2,
%%% for all t
%%% t = 1 is the initial period at time 0, final period t=101 is time 100
%%% we are at steady state at both 0 and 100, thus at times t=1 and t=101

% Exogenous monetary policy shock
shock_taylor = zeros([1, par.T+1]);

% Exogenous path for productivity
Z_path = ones([1, par.T+1]);

% Guess on the path for inflation
pi_guess = linspace(0, 0, par.T+1);

% Inflation shock for the Jacobian
shock_pi = zeros([1, par.T+1]);

% Initializing the Jacobian --> for time 1 to 100, thus t=2 to 101
Jacobian = zeros(par.T, par.T); 

% Initial distribution
init_dist = dist;


% Computing at steady state
simul = initialize(pi_guess, shock_pi, shock_taylor, Z_path, par);
Anext_matrix = backward_step(V, Anext, simul, par);
dist_matrix = forward_step(Anext_matrix, init_dist, par);
Asset_evolution = asset_evolution(dist_matrix, par);

path_steady_state = Asset_evolution;
% path_steady_state = ones(size(path_steady_state))*path_steady_state(end);




%% Computing the Jacobian 


for i = 2:par.T+1
    shock_pi = zeros([1, par.T+1]);
    shock_pi(i) = par.epsilon; 
    
    simul = initialize(pi_guess, shock_pi, shock_taylor, Z_path, par);
    Anext_matrix = backward_step(V, Anext, simul, par);
    dist_matrix = forward_step(Anext_matrix, init_dist, par);
    Asset_evolution = asset_evolution(dist_matrix, par);

    Jacobian(:, i-1) = (Asset_evolution - path_steady_state)' / par.epsilon;
    fprintf('\rJacobian Iteration %d', i); 

end 



%% Computing the transition to a one-time shock

% setting shock to inflation back to 0
shock_pi = zeros([1, par.T+1]);

% monetary policy shock
monetary_policy_shock = -0.0025; 
shock_taylor(1) = 1;
shock_taylor(2) = exp(monetary_policy_shock); 

for i =3:par.T+1
    shock_taylor(i) = exp(par.rho_nu * log(shock_taylor(i-1)));
end

shock_taylor = log(shock_taylor);

% persistence 
par.rho_nu = 0.5;

% Guess on the path for inflation
pi_guess = linspace(1e-2, 0, par.T+1) ;
pi_guess(1) = 0;

diff_metric = 1;
while diff_metric > 1e-12
    % solve with the updated guess
    simul = initialize(pi_guess, shock_pi, shock_taylor, Z_path, par);
    Anext_matrix = backward_step(V, Anext, simul, par);
    dist_matrix = forward_step(Anext_matrix, init_dist, par);
    Asset_evolution = asset_evolution(dist_matrix, par);
    
    % compute the residuals
    resid = Asset_evolution - path_steady_state;
    
    % update the guess
    pi_guess(2:par.T+1) = pi_guess(2:par.T+1) - (Jacobian \ resid')';
    % pi_guess = max(pi_guess, 0);
    
    % update the difference
    diff_metric = max(abs(Asset_evolution - path_steady_state));
    diff_metric
end



Consumption_path = (1-simul.tax_path) .* simul.wage_path .* simul.L_hh + (1+r)*B - [B, Asset_evolution]; 


%% Plotting the results 

xmin = 2;
xmax = 10;
figure;

% Plot Asset_evolution
subplot(3, 3, 1); % 3 rows, 2 columns, first plot
plot(Asset_evolution, 'b'); % 'b' specifies blue color
title('Aggregate Asset Demand');
xlabel('Time');
ylabel('Assets');
ylim([0,4])
xlim([xmin, xmax])
grid on;

% Plot pi_guess
subplot(3, 3, 2); % 3 rows, 2 columns, second plot
plot(pi_guess, 'b'); % 'b' specifies blue color
title('Inflation');
xlabel('Time');
ylabel('Inflation Rate');
xlim([xmin, xmax])
grid on;

% Plot r_path
subplot(3, 3, 3); % 3 rows, 2 columns, third plot
plot(simul.r_path, 'b'); % 'b' specifies blue color
title('Real Interest Rate');
xlabel('Time');
ylabel('Real Rate');
xlim([xmin, xmax])
grid on;

% Plot i_path
subplot(3, 3, 4); % 3 rows, 2 columns, fourth plot
plot(simul.i_path, 'b'); % 'b' specifies blue color
title('Monetary Policy');
xlabel('Time');
ylabel('Nominal Rate');
xlim([xmin, xmax])
grid on;

% Plot L_hh
subplot(3, 3, 5); % 3 rows, 2 columns, fifth plot
plot(simul.L_hh, 'b'); % 'b' specifies blue color
title('Labor Supply');
xlabel('Time');
ylabel('Labor');
xlim([xmin, xmax])
grid on;

% Plot tax_path
subplot(3, 3, 6); % 3 rows, 2 columns, sixth plot
plot(simul.tax_path, 'b'); % 'b' specifies blue color
title('Tax Level');
xlabel('Time');
ylabel('Tax');
xlim([xmin, xmax])
grid on;

% Plot consumption
subplot(3, 3, 7); % 3 rows, 2 columns, sixth plot
plot(Consumption_path, 'b'); % 'b' specifies blue color
title('Agg. Consumption');
xlabel('Time');
ylabel('Consumption');
xlim([xmin, xmax])
grid on;

% Plot consumption
subplot(3, 3, 8); % 3 rows, 2 columns, sixth plot
plot(shock_taylor, 'b'); % 'b' specifies blue color
title('Monetary P. Shock');
xlabel('Time');
ylabel('Shock');
xlim([xmin, xmax])
grid on;

% Adjust the layout to prevent overlap
sgtitle('Simulation Results');



%%
figure; 

it=0;
for i=1:9
    it = it+1;
    subplot(3, 3, it);
    
    di = dist_matrix(:,i);
    distribution = reshape(di, par.M, par.N);
    distribution = sum(distribution); 

    plot(par.Agrid, distribution, LineWidth=1, LineStyle="-", Color="blue")
    title(sprintf('Period %d', i));
    xlim([0, 10])
    ylim([0, 0.05])
    grid on


end

sgtitle('Distribution Evolution');