function K_Delta_tau = K_Delta_tau_func(l_over_lb, vr_over_cs, Param)
% K_DELTA_TAU_FUNC
% Background stress intensity factor
%
% Inputs:
%   l_over_lb        - crack length l / l_b
%   vr_over_cs       - rupture velocity v_r / c_s
%   bar_v0_over_cs   - scaled reference velocity \bar{v}_0
%   rs_type          - "slip" or "aging"
%
% Output:
%   K_Delta_tau      - background stress intensity factor (dimensionless)

    % Effective stress
    Delta_tau_eff = Delta_tau_eff_func( ...
        l_over_lb, vr_over_cs, Param);

    % Stress intensity factor
    K_Delta_tau = Delta_tau_eff .* sqrt(pi .* l_over_lb);
end
