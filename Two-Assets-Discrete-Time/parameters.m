function par  = parameters()

%% Parameters

% Households
par.gamma    = 2;
par.eis = 1/par.gamma;
par.beta  = 0.98;

% Labour
par.M       = 3;
par.rho     = 0.95;
par.sigma_eps   = 0.05;

% Grids size
par.nba      = 50;
par.nbb      = 50;

% Firm
par.Z       = 1;

% Government
par.debt = 2;

% NKPC
par.mu = 2;
par.kappa = 0.1;
par.phi = 1.5;
par.frisch = 1;

% Illiquid Asset grid
par.Amin   = 0;
par.Amax   = 100;
ubar = log(1 +log(1 + par.Amax - par.Amin));
u_grid = linspace(0, ubar, par.nba);
par.Agrid = par.Amin + exp(exp(u_grid) - 1) - 1;
par.Agrid = linspace(par.Amin, par.Amax, par.nba);


% Liquid Asset grid
par.Bmin   = 0;
par.Bmax   = 100;
ubar = log(1 +log(1 + par.Bmax - par.Bmin));
u_grid = linspace(0, ubar, par.nbb);
par.Bgrid = par.Bmin + exp(exp(u_grid) - 1) - 1;
par.Bgrid = linspace(par.Bmin, par.Bmax, par.nbb);

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
par.w = par.labor;

% Big grids

[par.SS,par.AA] = ndgrid(par.s,par.Agrid);


[par.AAgrid, par.BBgrid] = meshgrid(par.Agrid, par.Bgrid);
[par.BBB, par.AAA, par.SSS] = ndgrid(par.Bgrid, par.Agrid, par.s);



par.S = [par.AA(:),par.SS(:)];
par.PP = kron(speye(par.nba),par.P);


% Adjustment cost function 

par.g = @(a, ap) 0.1 * (a ~= ap);
par.g = @(a, ap)  (a - ap).^2;


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