function [Anext, Bnext, V] = vfi_Graves(par)

% Parameters
eps = 1e-6;        
diff_metric = 1;          
n_step = 0;         
max_iter = 1000;    


V = (par.AAA + par.BBB + par.SSS).^(par.eis);

% Preallocate for speed
V_int = zeros(par.nbb, par.nba, par.M);
A_pol = zeros(par.nbb, par.nba, par.M);
Ev = zeros(par.nbb, par.nba, par.M);




%% First step 

while diff_metric > eps && n_step < max_iter
    n_step = n_step + 1;
    
    % Reshape V to a 2D matrix: (nbb*nba) x M
    V_reshaped = reshape(V, [], par.M);  % size: (nbb*nba) x M
    % Compute expected value: weighted sum over future states
    Ev_reshaped = par.beta * (V_reshaped * par.P');  % size: (nbb*nba) x M
    % Reshape back to original shape
    Ev = reshape(Ev_reshaped, par.nbb, par.nba, par.M);


    V_interp = griddedInterpolant(par.BBB, par.AAA, par.SSS, Ev,'linear', 'linear');
    
    [B_pol_NA, V_int] = golden_search_NA(V_interp, par);

   
    % Checking for convergence
    diff_metric = max(max(max(abs(V_int-V)))); 

    % Display the location
    % abs_diff = abs(V - V_int);
    % max_value = max(abs_diff(:));
    % [iz, ib, ik] = ind2sub(size(abs_diff), find(abs_diff == max_value, 1));
    % fprintf('The location of the maximum absolute difference is at (z, b, k) = (%d, %d, %d)\n', iz, ib, ik);

    if diff_metric < eps
        disp(['First step converged in ', num2str(n_step), ' iterations'])
    end
    disp(diff_metric)
    
    % Update the value function
    V = V_int; 
end 

par.Bnext_NA = B_pol_NA;
par.Bnext_NA_grid = griddedInterpolant(par.BBB, par.AAA, par.SSS, B_pol_NA, 'linear', 'linear');

temp = griddedInterpolant(par.BBB, par.AAA, par.SSS, B_pol_NA, 'linear', 'linear');



%% Second step


diff_metric = 1; 
n_step = 0;
V = (par.AAA + par.BBB + par.SSS).^(par.eis);
max_iter = 1000;

while diff_metric > eps && n_step < max_iter
    n_step = n_step + 1;
    
    % Reshape V to a 2D matrix: (nbb*nba) x M
    V_reshaped = reshape(V, [], par.M);  % size: (nbb*nba) x M
    % Compute expected value: weighted sum over future states
    Ev_reshaped = par.beta * (V_reshaped * par.P');  % size: (nbb*nba) x M
    % Reshape back to original shape
    Ev = reshape(Ev_reshaped, par.nbb, par.nba, par.M);


    V_interp = griddedInterpolant(par.BBB, par.AAA, par.SSS, Ev,'linear', 'linear');
    
    [A_pol, B_pol, V_int] = golden_search_Adj(V_interp, par);

   
    % Checking for convergence
    diff_metric = max(max(max(abs(V_int-V)))); 

    % Display the location
    % abs_diff = abs(V - V_int);
    % max_value = max(abs_diff(:));
    % [iz, ib, ik] = ind2sub(size(abs_diff), find(abs_diff == max_value, 1));
    % fprintf('The location of the maximum absolute difference is at (z, b, k) = (%d, %d, %d)\n', iz, ib, ik);

    if diff_metric < eps
        disp(['Second step converged in ', num2str(n_step), ' iterations'])
    end
    disp(diff_metric)
    
    % Update the value function
    V = V_int; 
end 

Anext = A_pol;
Bnext = B_pol;



end