function [f g J] = misfit(m, Do, Q, model);

    [Dp, J] = F(m,Q,model);
    dD      = Dp - Do;
    g       = J'*dD;
    dD      = gather(dD);
    f       = 0.5 * norm(dD)^2;
