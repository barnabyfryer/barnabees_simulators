function Veff = Veff_func(l_over_lb, vr_over_cs, Param)
% VEFF_FUNC
% Effective slip velocity entering the rate-and-state law
%
% Inputs:
%   l_over_lb        - crack length l / l_b
%   vr_over_cs       - rupture velocity v_r / c_s
%   bar_v0_over_cs   - scaled reference velocity \bar{v}_0
%   rs_type          - "slip" or "aging"
%
% Output:
%   Veff             - effective velocity (dimensionless)

    % --- Constants
    F_cal = 1;

    % --- Elastic functions
    g = g_func(vr_over_cs);
    k = k_func(vr_over_cs);

    % --- Scaled toughness
    Kc = Kc_func(vr_over_cs, Param);

    % --- Effective velocity
    Veff = (k ./ g) ...
        .* (4 .* Kc .* vr_over_cs ./ Param.v0_over_cs) ...
        ./ sqrt(pi .* l_over_lb) ...
        .* F_cal;
end
