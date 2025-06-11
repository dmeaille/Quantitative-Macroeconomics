function [A, dist, anext, C, V] = household_ss(r, V, par)

[anext, C, V] = egm(r, V, par);
T = compute_Trans_matrix_transp(anext, par); 


% Initial distribution
dist=zeros(par.nba*par.nb_states,1);
% Initialize mass at the middle
dist(ceil(par.nba*par.nb_states/2))=1;
error=1;
it=1;
while error>1e-8
    dist1=T'*dist;
    error=max(max(abs(dist1-dist)));
    dist=dist1;
    it=it+1;
end

A=dist'*kron(par.agrid',ones(3,1));

end

