function [distribution_global, distribution_by_states, c_policy] = general_equilibrium(par)

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


distribution_by_states = reshape(dist, par.nb_states, par.nba);
distribution_global = sum(distribution_by_states);


end