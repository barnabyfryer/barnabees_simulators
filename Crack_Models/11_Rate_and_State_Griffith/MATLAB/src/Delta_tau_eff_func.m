function Delta_tau_eff = Delta_tau_eff_func( ...
        l_over_lb, vr_over_cs, Param)
% DELTA_TAU_EFF_FUNC
% Effective stress drop entering the EoM
%
% Inputs:
%   l_over_lb        - crack length l / l_b
%   vr_over_cs       - rupture velocity v_r / c_s
%   bar_v0_over_cs   - scaled reference velocity \bar{v}_0
%   rs_type          - "slip" or "aging"
%
% Output:
%   Delta_tau_eff    - effective stress drop

    % --- Effective velocity ratio
    Veff_over_V0 = Veff_func(l_over_lb, vr_over_cs, Param);

    % --- Effective stress
    Delta_tau_eff = Param.Delta_f0_over_b ...
        - (Param.a_over_b - 1) .* log(Veff_over_V0);
end
