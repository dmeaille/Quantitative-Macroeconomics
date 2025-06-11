function Vnew = vfun_nested(Bprime, Aprime, par, V_interp)
    % b in row, a in column
    % initialization 
    c = (1+par.r_a)*par.AAA + (1+par.r_b) * par.BBB - par.g(par.AAA, Aprime) + par.SSS*par.w - Aprime - Bprime;
    u = ones(par.nbb, par.nba, par.M)*-1e5; 
    Vnew = ones(size(Bprime))*-1e5; 

    % excluding c<0
    u = -1e5 * ones(size(c)); 
    valid = c > 0;
    u(valid) = (c(valid).^(1 - par.gamma)) / (1 - par.gamma);
    
    % computing values for next period value function 
    
    % need to interp also for kprime interp1(param.bgrid, Ev(i,:), bprime(i,:), 'spline', 'extrap');
    % computing U + B*V' (the beta discounting factor is factored in the
    % interpolator !
    Vnew = u + V_interp(Bprime, Aprime, par.SSS);
    
end