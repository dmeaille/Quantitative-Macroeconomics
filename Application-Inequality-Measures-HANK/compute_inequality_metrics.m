function [gini, asset_ratio, consumption_ratio] = compute_inequality_metrics(par)


[dist_global, dist_states, c_policy] = general_equilibrium(par);



% Gini Index

cdf = cumsum(dist_global);

asset_total = sum(par.agrid .* dist_global);
asset_share_cum = cumsum(par.agrid .* dist_global) / asset_total;

cdf_lorenz = [0, cdf];
asset_lorenz = [0, asset_share_cum]; % Lorenz curve 

% Gini is 1 - 2 * area under Lorenz curve
gini = 1 - 2 * trapz(cdf_lorenz, asset_lorenz);


% Asset Ratio (par.agrid is already in increasing order)

% Top 10% vs Bottom 10%
top10_idx = (cdf >= 0.9);  % logical index of top 10%
bottom10_idx = (cdf <= 0.1); % logical index of bottom 10%

asset_top10 = sum(par.agrid(top10_idx) .* dist_global(top10_idx));
asset_bottom10 = sum(par.agrid(bottom10_idx) .* dist_global(bottom10_idx));

asset_ratio = asset_top10 / asset_bottom10;



% Consumption Ratio 

% Ordering the consumption choices with the distribution allocated to each
c_policy = c_policy(:);
dist_states = dist_states(:);

[c_pol_sorted, indx] = sort(c_policy);
da_sorted = dist_states(indx);

cdf_c = cumsum(da_sorted);

% Top 10% vs Bottom 10%
top10_idx = (cdf_c >= 0.9);  % logical index of top 10%
bottom10_idx = (cdf_c <= 0.1); % logical index of bottom 10%

cons_top10 = sum(c_pol_sorted(top10_idx) .* da_sorted(top10_idx));
cons_bottom10 = sum(c_pol_sorted(bottom10_idx) .* da_sorted(bottom10_idx));

consumption_ratio = cons_top10 / cons_bottom10;


end