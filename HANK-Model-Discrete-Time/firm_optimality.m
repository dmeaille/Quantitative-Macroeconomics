function [pi_wage_path] = firm_optimality(pi_path, Z_path, par)

    
    pi_wage_path = pi_path;
    
    % pi_wage_path(1) = pi_path(1);
    % pi_wage_path(2:par.T-1) = (1+pi_path(2:par.T-1)) .* Z_path(2:par.T-1)./Z_path(1:par.T-2) - 1;
    % pi_wage_path(end) = 0;
    

end