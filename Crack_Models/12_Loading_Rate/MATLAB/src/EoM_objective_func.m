function res = EoM_objective_func(l_over_lb, vr_over_cs, Param)
% EOM_OBJECTIVE_FUNC
% Residual of the equation of motion for rupture evolution.
%
% Inputs:
%   l_over_lb       - crack length (l / l_b)
%   vr_over_cs      - rupture velocity (v_r / c_s)
%   bar_v0_over_cs  - scaled ambient rupture velocity
%   rs_type         - 'slip' or 'aging' friction law
%
% Output:
%   res             - residual value: K_Delta_tau + K_Delta_T - Kc

    % Fracture toughness
    Kc = Kc_func(vr_over_cs, Param);
    
    % Background stress contribution
    K_Delta_tau = K_Delta_tau_func(l_over_lb, vr_over_cs, Param);
    
    % Residual of EoM
    res = K_Delta_tau - Kc;

end
