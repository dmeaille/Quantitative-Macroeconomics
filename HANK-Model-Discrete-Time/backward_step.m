function [Anext_matrix] = backward_step(V, Anext, simul, par)

% Inverted EGM, find previous period as a function of today
Vnext = V; 
Anext_matrix = zeros(par.T+1, par.M, par.N);
Anext_matrix(par.T+1,:,:) = Anext; 

for t=par.T:-1:1
    [A_prime, ~, V_temp] = egm_step(simul.r_path(t), simul.wage_path(t), simul.tax_path(t), Vnext, simul.L_hh(t), par);
    Anext_matrix(t,:,:) = A_prime;
    Vnext = V_temp;
end
