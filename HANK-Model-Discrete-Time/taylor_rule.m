function [i_path] = taylor_rule(pi_path, shock_taylor, par)

    i_path = par.r + 1.5 * pi_path  + shock_taylor;
    %i_path = (1 + par.r) * pi_path.^1.5  .* shock_taylor;

end