function [Anext, C, V] = egm(r, w, tax, V, N, par)

error = 1;
Aold = par.AA;


while error > par.tol_egm

    RHS = par.beta * par.P * V;

    C = RHS.^(-par.eis);
    coh = par.AA * (1+r) + (1-tax) * w * par.SS * N; % cash on hand, N is indiv labor supply

    Anext = zeros(par.M, par.N);

    for is = 1:par.M
        Anext(is,:) = interp1(par.Agrid + C(is,:), par.Agrid, coh(is,:),"linear","extrap");
        % Here: write your own & fast vectorized linear interpolation or
        % use griddedInterpolant
    end
    %grid = griddedInterpolant(Aold + C, Aold, 'linear', 'extrap'); % Define the interpolant
    %Anext(is, :) = grid(coh(is, :));

    Anext(Anext < par.Agrid(1)) = par.Agrid(1);
    Anext(Anext > par.Agrid(end)) = par.Agrid(end) - 1e-4;
    C = coh - Anext;

    V = (1+r) * C.^(-1/par.eis);
   
    % [Anext, C, V] = egm_step(r, w, tax, 1, V, par);
    error = sum(abs(Anext(:) - Aold(:)));
    Aold = Anext;
end


end


