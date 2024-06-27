function [residual, g1, g2] = dynamic_resid_g1_g2(T, y, x, params, steady_state, it_, T_flag)
% function [residual, g1, g2] = dynamic_resid_g1_g2(T, y, x, params, steady_state, it_, T_flag)
%
% Wrapper function automatically created by Dynare
%

    if T_flag
        T = DCP_D.dynamic_g2_tt(T, y, x, params, steady_state, it_);
    end
    [residual, g1] = DCP_D.dynamic_resid_g1(T, y, x, params, steady_state, it_, false);
    g2       = DCP_D.dynamic_g2(T, y, x, params, steady_state, it_, false);

end
