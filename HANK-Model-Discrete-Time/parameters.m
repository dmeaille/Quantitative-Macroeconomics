function par  = parameters()

%% Parameters

% Households
par.eis    = 0.5;
par.gamma = 1/par.eis;
par.beta  = 0.95;

% Labour
par.M       = 5;
par.rho     = 0.9;
par.sigma_eps   = 0.2;

% Firm
par.Z       = 1;

% Government
par.debt = 2;

% NKPC
par.mu = 2;
par.kappa = 0.1;
par.phi = 1.5;
par.frisch = 1;

% Asset grid
par.N      = 250;
par.Amin   = 0;
par.Amax   = 300;

ubar = log(1 +log(1 + par.Amax - par.Amin));
u_grid = linspace(0, ubar, par.N);
par.Agrid = par.Amin + exp(exp(u_grid) - 1) - 1;

% Invariant distribution of labor

[logs,par.P] = tauchen(par.M,0,par.rho,par.sigma_eps,1.8); 
logs         = logs';


invdist = ones(1,par.M)/par.M; 
test    = 1; 

while (test>0.0000001)
    invdist1 = invdist*par.P;
    test     = max(abs(invdist1-invdist));
    invdist  = invdist1;
end

par.invdist     = invdist';
par.s           = exp(logs);

labor = par.s * par.invdist;
par.s = par.s / labor;
par.labor = par.s * par.invdist;

% Big grids

[par.AA,par.SS] = meshgrid(par.Agrid,par.s);
par.S = [par.AA(:),par.SS(:)];
par.PP = kron(speye(par.N),par.P);

% Algorithm

par.tol_clear = 1e-8;
par.tol_egm   = 1e-8;
par.tol_trans = 1e-8;
par.tol_bis = 1e-10;
par.r_guess = 0.02; 

par.epsilon = 1e-5; % For the derivative

% Transition
par.T = 100;
% Shock
par.rho_nu = 0.5;
end