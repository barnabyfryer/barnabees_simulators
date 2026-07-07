function K_Delta_T = K_Delta_T_func(l_over_lb, Delta_T)
% K_DELTA_T_FUNC
% Foreshock (hypocentral) stress intensity factor
%
% Inputs:
%   l_over_lb   - crack length l / l_b
%   Delta_T    - scaled hypocentral force
%
% Output:
%   K_Delta_T  - foreshock stress intensity factor (dimensionless)

    K_Delta_T = Delta_T ./ sqrt(pi .* l_over_lb);
end
