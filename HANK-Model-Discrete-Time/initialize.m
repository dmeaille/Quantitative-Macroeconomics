function [simul] = initialize(pi_guess, shock_pi, shock_taylor, Z_path, par)


pi_path = pi_guess + shock_pi; 

% Induced paths

i_path = taylor_rule(pi_path, shock_taylor, par);
r_path = fisher(i_path, pi_path, par);
%r_path = ones([1, par.T])*par.r;

pi_wage_path = firm_optimality(pi_path, Z_path, par); 
wage_path = zeros([1, par.T+1]); 
wage_path = Z_path;
% wage_path(1) = 1+(pi_wage_path(1))*par.Z;
% for t=2:par.T+1
%     wage_path(t) = (1+pi_wage_path(t))*wage_path(t-1);
% end


% Finding Aggregate Labour Supply
L0 = ones([1, par.T]);
L_hh = nkpc_find_l(pi_wage_path, wage_path, r_path, L0, par); 

% Implied Tax rates
tax_path(1:par.T) = (par.debt*r_path(1:par.T))./(wage_path(1:par.T).*L_hh(1:par.T));
tax_path(par.T+1) = par.tax; 



simul.i_path = i_path; 
simul.r_path = r_path; 
simul.pi_wage_path = pi_wage_path;
simul.wage_path = wage_path; 
simul.L_hh = L_hh;
simul.tax_path = tax_path; 






end