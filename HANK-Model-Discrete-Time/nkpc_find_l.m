function [L] = nkpc_find_l(pi_wage_path, wage_path, r_path, L0, par)

h = 1e-5;
f_x = 1;
nn = size(L0);

if nn(2) == par.T 
    L = L0;
else
    L = L0(1:par.T);
end

nmin = 1/3;
nmax = 1;
w = wage_path(1:par.T);
piw_t1 = pi_wage_path(2:par.T+1);
piw_t = pi_wage_path(1:par.T);

while max(abs(f_x)) > 0.0001

    tax = (par.debt*r_path(1:par.T))./(L);
    tax_bis = (par.debt*r_path(1:par.T))./(L+h);
    
    %- 1/par.mu * (1-par.r*par.debt/nmin) * w * nmin^(-1/par.eis) ) ...
    
    
    f_x = par.kappa * ( par.phi * L.^(par.frisch) ...
        - 1/par.mu * (1 - tax) .* w .* (L.*w).^(-par.gamma)) ...
        + par.beta * piw_t1  - piw_t;

    %- 1/par.mu * (1-par.r*par.debt/nmin) * w * nmin^(-1/par.eis)) ...

    fprime_x = (par.kappa * ( par.phi * (L+h).^(par.frisch) ...
        - 1/par.mu * (1 - tax_bis) .* w .* ((L+h).*w).^(-par.gamma)) ...
        + par.beta * piw_t1  - piw_t - f_x)./h;
    
    L_new = L - f_x./fprime_x;
    L = L_new;

end

L(par.T+1) = 1;
% L = min(L, nmax);
% L = max(L, nmin);

end