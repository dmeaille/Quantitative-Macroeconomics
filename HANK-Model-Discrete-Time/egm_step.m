function [Anext, C, V] = egm_step(r, w, tax, V, N, par)
   
    RHS = par.beta * par.P * V;                
    
    C = RHS.^(-par.eis);
    coh = par.AA * (1+r) + (1-tax) * w * par.SS * N;
    
    for is = 1:par.M
        Anext(is,:) = interp1(par.Agrid + C(is,:), par.Agrid, coh(is,:),"linear","extrap");
        % Here: write your own & fast vectorized linear interpolation or
        % use griddedInterpolant
    end
            
    Anext(Anext < par.Agrid(1)) = par.Agrid(1);
    Anext(Anext > par.Agrid(end)) = par.Agrid(end) - 1e-4;
    C = coh - Anext;
    
    V = (1+r) * C.^(-1/par.eis);
end
        