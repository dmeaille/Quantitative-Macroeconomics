function [dist_matrix] = forward_step(Anext_matrix, init_dist, par)

dist_matrix = zeros(par.M*par.N, par.T+1); 
dist_matrix(:, 1) = init_dist; 

for t=2:par.T+1
    A_prime = squeeze(Anext_matrix(t, :, :));
    Trans = compute_Trans_matrix_transp(A_prime, par);
    dist_matrix(:, t) = Trans' * dist_matrix(:, t-1); 
end