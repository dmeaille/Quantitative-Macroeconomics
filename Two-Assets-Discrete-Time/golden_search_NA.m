function [b_p, f_bp] = golden_search_NA(V_interp, par)


% Initialization
psy = (sqrt(5) + 1)/2; 
alpha = 1/psy; 

max_iter = 1000; 
count_iter = 0; 
crit = 1e-7; 
next=1;

% Algorithm
h_u = min(par.BBB+1, par.Bmax); 
h_l = zeros(size(par.BBB));


for i=1:max_iter

    % Compute test points
    b_l = h_l + (1 - alpha) .* (h_u - h_l);
    b_u = h_u - (1 - alpha) .* (h_u - h_l);

    % Evaluate function at test points --> no adjustment on the capital
    % decision
    fun_u = vfun_TwoAssets(b_u, par.AAA, par, V_interp); 
    fun_l = vfun_TwoAssets(b_l, par.AAA, par, V_interp); 

    % Update bounds using logical indexing
    mask = fun_u > fun_l;
    h_l(mask) = b_l(mask);
    h_u(~mask) = b_u(~mask);

    % Convergence check for entire grid
    if max(max(abs(fun_u-fun_l))) < crit %all(abs(fun_u - fun_l) < crit, 'all')
        break;
    end

    count_iter = count_iter + 1; 
end

b_p = (h_u + h_l)/2;
f_bp = vfun_TwoAssets(b_p, par.AAA, par, V_interp);


end