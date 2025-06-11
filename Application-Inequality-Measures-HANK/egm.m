function [anext, C, V] = egm(r, V, par)

error = 1;
Aold = par.Amesh;


while error > par.tol_egm

    RHS = par.beta * par.P * V;

    C = RHS.^(-par.eis);
    coh = par.Amesh * (1+r) + par.Ymesh; % cash on hand

    anext = zeros(par.nb_states, par.nba);

    for is = 1:par.nb_states
        anext(is,:) = interp1(par.agrid + C(is,:), par.agrid, coh(is,:),"linear","extrap");
    end


    anext(anext < par.agrid(1)) = par.agrid(1);
    anext(anext > par.agrid(end)) = par.agrid(end) - 1e-4;
    C = coh - anext;

    V = (1+r) * C.^(-1/par.eis);
   
    % [anext, C, V] = egm_step(r, w, tax, 1, V, par);
    error = sum(abs(anext(:) - Aold(:)));
    Aold = anext;
end


end


