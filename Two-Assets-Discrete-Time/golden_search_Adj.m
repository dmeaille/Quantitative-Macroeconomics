function [a_p, b_p, f_p] = golden_search_Adj(V_interp, par)


% Initialization
psy = (sqrt(5) + 1)/2; 
alpha = 1/psy; 

max_iter = 1000; 
count_iter = 0; 
crit = 1e-8; 
next=1;

% Algorithm
h_u = min(par.AAA+1, par.Amax+1); 
h_l = zeros(size(par.AAA));


for i=1:max_iter

    % Compute test points
    a_l = h_l + (1 - alpha) .* (h_u - h_l);
    a_u = h_u - (1 - alpha) .* (h_u - h_l);

    % Find implied bond adjustment policy
    b_star_l = par.BBB + (1+par.r_a)/(1+par.r_b) * (par.AAA - a_l) - par.g(par.AAA, a_l)/(1+par.r_b);
    b_star_u = par.BBB + (1+par.r_a)/(1+par.r_b) * (par.AAA - a_u) - par.g(par.AAA, a_u)/(1+par.r_b);


    b_adj_l = par.Bnext_NA_grid(b_star_l, a_l, par.SSS);
    b_adj_u = par.Bnext_NA_grid(b_star_u, a_u, par.SSS);

    % b_adj_l = zeros(par.nbb, par.nba, par.M);
    % b_adj_u = zeros(par.nbb, par.nba, par.M);
    % M = par.M;
    % parfor i =1:M
    %     b_adj_l(:,:,i) = interp2(par.AAgrid, par.BBgrid, squeeze(par.Bnext_NA(:,:,i)), squeeze(a_l(:,:,i)), squeeze(b_star_l(:,:,i)), 'spline');
    %     b_adj_u(:,:,i) = interp2(par.AAgrid, par.BBgrid, squeeze(par.Bnext_NA(:,:,i)), squeeze(a_u(:,:,i)), squeeze(b_star_u(:,:,i)), 'spline');
    % end

    % Evaluate function at test points
    fun_u = vfun_TwoAssets(b_adj_u, a_u, par, V_interp); 
    fun_l = vfun_TwoAssets(b_adj_l, a_l, par, V_interp); 

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

% when adjusting
a_p = (h_u + h_l)/2; 
b_p_star = par.BBB + (1+par.r_a)/(1+par.r_b) * (par.AAA - a_l) - par.g(par.AAA, a_l)/(1+par.r_b);
b_p = par.Bnext_NA_grid(b_p_star, a_p, par.SSS);
f_p = vfun_TwoAssets(b_p, a_p, par, V_interp);


if next
    % when not adjusting
    b_p_no_adj = par.Bnext_NA;
    f_no_adj = vfun_TwoAssets(b_p_no_adj, par.AAA, par, V_interp);

    % mask to know when we adjust
    mask = f_no_adj > f_p; 
    
    % Gumbel shock 
    % f_max = max(f_no_adj, f_bp); 
    % exp_sum = exp((f_no_adj - f_max)./par.gumbel_par) + exp((f_bp - f_max)./par.gumbel_par); 
    % 
    % prob_adj = exp((f_bp - f_max)./par.gumbel_par)./exp_sum; 
    % prob_no_adj = exp((f_no_adj - f_max)./par.gumbel_par)./exp_sum; 

    % updating the outputs
    a_p(mask) = par.AAA(mask);
    b_p(mask) = b_p_no_adj(mask);
    f_p = vfun_TwoAssets(b_p, a_p, par, V_interp); 
    
    % mask means: we adjust if 1
    mask = ~mask; 
end

end