function [Asset_evolution] = asset_evolution(dist_matrix, par)

Asset_evolution = zeros(1, par.T); 

for t=1:par.T 
    distt = squeeze(dist_matrix(:, t+1));
    Asset_evolution(t) = distt'*kron(par.Agrid',ones(5,1));
end

