function [r_path] = fisher(i_path, pi_path, par)


    r_path = (1+i_path(1:par.T))./(1+pi_path(2:par.T+1)) - 1;
    r_path(par.T+1) = i_path(par.T+1);
    
    %r_path(par.T) = par.r; 
    %r_path = (1+i_path(1:par.T))./(1+pi_path(2:end)) - 1;


end