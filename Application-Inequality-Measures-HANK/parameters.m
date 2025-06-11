function par  = parameters()

%% Parameters

% Households
par.eis    = 0.5;
par.gamma = 1/par.eis;
par.beta  = 0.95;

% Labour
par.nb_states       = 3;

% Firm TFP
par.Z       = 1;

% Firm Asset
par.capital = 2;
par.alpha = 0.34;


% Asset grid
par.nba      = 500;
par.Amin   = 0;
par.Amax   = 20;

ubar = log(1 +log(1 + par.Amax - par.Amin));
u_grid = linspace(0, ubar, par.nba);
par.agrid = par.Amin + exp(exp(u_grid) - 1) - 1;

% Invariant distribution of employment

pi_g_g = 0.9;
pi_b_b = 0.8; 
pi_u_g = 0.3; 
pi_u_b = 0.5; 

% transition matrix
P = [[pi_g_g, 0, 1-pi_g_g];
    [0, pi_b_b, 1-pi_b_b];
    [pi_u_g, pi_u_b, 1-pi_u_b-pi_u_g]];

par.P = P;


% steady-state distribution 
invdist = ones(1,par.nb_states)/par.nb_states; 
test    = 1; 

while (test>0.0000001)
    invdist1 = invdist*par.P;
    test     = max(abs(invdist1-invdist));
    invdist  = invdist1;
end

par.invdist     = invdist';



tau = 0.3;      % tax rate 
mu_g = 0.95;
mu_b = 0.85;
mu_u = 0.3;

par.Gamma = mu_g/mu_b-1; 
par.tau = tau; 
par.mu_g = (1+par.Gamma)*mu_b;
par.mu_b = mu_b;
par.mu_u = mu_u;

% vector of mus
mu = [mu_g, mu_b, mu_u];

N_g = invdist(1);
N_b = invdist(2);

N = [N_g, N_b, 0];

par.labor = mu * N';
% labor = par.y .* par.invdist;
% par.y = par.y ./ labor;
% par.labor = par.y * par.invdist;




% Wage level
w = 1;
par.wage=w;
income = [par.mu_g * par.wage * (1-par.tau), par.mu_b * par.wage * (1-par.tau), par.mu_u * par.wage];
% Z = Z';     % income 

par.y = income;


% Big grids

[par.Amesh,par.Ymesh] = meshgrid(par.agrid,par.y);
% par.S = [par.Amesh(:),par.Ymesh(:)];
% par.PP = kron(speye(par.nba),par.P);

% Algorithm

par.tol_clear = 1e-8;
par.tol_egm   = 1e-8;
par.tol_trans = 1e-8;
par.tol_bis = 1e-10;




end