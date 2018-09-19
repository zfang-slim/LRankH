function [f dD J] = misfit_GN(m, Do, Q, model);

    [Dp, J] = F(m,Q,model);
    dD      = Dp - Do;
    f       = 0.5 * gather(norm(dD)^2);
