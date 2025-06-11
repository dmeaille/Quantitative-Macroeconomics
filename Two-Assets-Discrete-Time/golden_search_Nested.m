function [b_p, f_bp] = golden_nested_find_bp(Aprime, par, V_interp)


% Initialization
psy = (sqrt(5) + 1)/2; 
alpha = 1/psy; 

max_iter = 30; 
count_iter = 0; 
crit = 1e-8; 
next=1;

% Algorithm
h_u = min(par.BBB + 1, par.Bmax); 
h_l = zeros(size(par.BBB));


for i=1:max_iter

    % Compute test points
    b_l = h_l + (1 - alpha) .* (h_u - h_l);
    b_u = h_u - (1 - alpha) .* (h_u - h_l);

    % Evaluate function at test points
    fun_u = vfun_TwoAssets(b_u, Aprime, par, V_interp); 
    fun_l = vfun_TwoAssets(b_l, Aprime, par, V_interp); 

    % Update bounds using logical indexing
    mask = fun_u > fun_l;
    h_l(mask) = b_l(mask);
    h_u(~mask) = b_u(~mask);

    % Convergence check for entire grid
    if max(max(abs(fun_u - fun_l))) < crit %all(abs(fun_u - fun_l) < crit, 'all')
        break;
    end

    count_iter = count_iter + 1; 
end

if next
    b_p = (h_u + h_l)/2; 
    f_bp = vfun_TwoAssets(b_p, Aprime, par, V_interp);
end