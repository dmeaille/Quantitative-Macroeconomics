function [a_p, b_p, f_bp] = golden_search_two_assets_nested(V_interp, par)


% Initialization
psy = (sqrt(5) + 1)/2; 
alpha = 1/psy; 

max_iter = 30; 
count_iter = 0; 
crit = 1e-8; 
next=1;

% Algorithm
h_u = min(par.AAA + 1, par.Amax+2); 
h_l = zeros(size(par.AAA));


for i=1:max_iter

    % Compute test points
    a_l = h_l + (1 - alpha) .* (h_u - h_l);
    a_u = h_u - (1 - alpha) .* (h_u - h_l);

    % Find corresponding optimal b policy
    b_policy_l = golden_search_Nested(a_l, par, V_interp);
    b_policy_u = golden_search_Nested(a_u, par, V_interp);

    % Evaluate function at test points
    fun_u = vfun_TwoAssets(b_policy_u, a_u, par, V_interp); 
    fun_l = vfun_TwoAssets(b_policy_l, a_l, par, V_interp); 

    % Update bounds using logical indexing
    mask = fun_u > fun_l;
    h_l(mask) = a_l(mask);
    h_u(~mask) = a_u(~mask);

    % Convergence check for entire grid
    if max(max(abs(fun_u-fun_l))) < crit %all(abs(fun_u - fun_l) < crit, 'all')
        break;
    end

    count_iter = count_iter + 1; 
end

if next
    % when adjusting
    a_p = (h_u + h_l)/2; 
    b_p = golden_search_Nested(a_p, par, V_interp);
    f_bp = vfun_TwoAssets(b_p, a_p, par, V_interp);

    % when not adjusting
    b_p_no_adj = golden_search_Nested(par.AAA, par, V_interp);
    f_no_adj = vfun_TwoAssets(b_p_no_adj, par.AAA, par, V_interp);

    % mask to know when we adjust
    mask = f_no_adj > f_bp; 
    
    % Gumbel shock 
    % f_max = max(f_no_adj, f_bp); 
    % exp_sum = exp((f_no_adj - f_max)./par.gumbel_par) + exp((f_bp - f_max)./par.gumbel_par); 
    % 
    % prob_adj = exp((f_bp - f_max)./par.gumbel_par)./exp_sum; 
    % prob_no_adj = exp((f_no_adj - f_max)./par.gumbel_par)./exp_sum; 

    % updating the outputs
    a_p(mask) = par.AAA(mask);
    b_p(mask) = b_p_no_adj(mask);
    f_bp = vfun_TwoAssets(b_p, a_p, par, V_interp); 
    
    % mask means: we adjust if 1
    mask = ~mask;

    % k_p = kprime .* prob_no_adj + k_p * prob_adj; 
    % k_b = b_p_no_adj .* prob_no_adj + b_p * prob_adj; 
    % f_bp = f_no_adj .* prob_no_adj + f_bp * prob_adj; 
end