function [Anext, Bnext, V] = vfi_TA(par)

% Parameters
eps = 1e-4;        
diff_metric = 1;          
n_step = 0;         
max_iter = 1000;    


V = (par.AAA + par.BBB + par.SSS).^(par.eis);

% Preallocate for speed
V_int = zeros(par.nbb, par.nba, par.M);
A_pol = zeros(par.nbb, par.nba, par.M);
Ev = zeros(par.nbb, par.nba, par.M);
tic
while diff_metric > eps && n_step < max_iter
    n_step = n_step + 1;
    
    % Reshape V to a 2D matrix: (nbb*nba) x M
    V_reshaped = reshape(V, [], par.M);  % size: (nbb*nba) x M
    % Compute expected value: weighted sum over future states
    Ev_reshaped = par.beta * (V_reshaped * par.P');  % size: (nbb*nba) x M
    % Reshape back to original shape
    Ev = reshape(Ev_reshaped, par.nbb, par.nba, par.M);


    V_interp = griddedInterpolant(par.BBB, par.AAA, par.SSS, Ev,'linear', 'linear');
    [A_pol, B_pol, V_int] = golden_search_TwoAssets(V_interp, par);

   
    % Checking for convergence
    diff_metric = max(max(max(abs(V_int-V)))); 

    % Display the location
    % abs_diff = abs(V - V_int);
    % max_value = max(abs_diff(:));
    % [iz, ib, ik] = ind2sub(size(abs_diff), find(abs_diff == max_value, 1));
    % fprintf('The location of the maximum absolute difference is at (z, b, k) = (%d, %d, %d)\n', iz, ib, ik);

    if diff_metric < eps
        disp(['Converged in ', num2str(n_step), ' iterations'])
    end
    disp(diff_metric)
    
    % Update the value function
    V = V_int; 
end 
toc

Anext = A_pol;
Bnext = B_pol;


end