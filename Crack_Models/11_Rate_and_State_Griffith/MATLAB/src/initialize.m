function vr_over_cs = initialize(l_over_lb, Param)
% INITIALIZE
% Find the initial rupture velocity v_r/c_s for a given crack length l/l_b
%
% Inputs:
%   l_over_lb      - initial crack length (l / l_b)
%   Delta_T        - scaled hypocentral force
%   bar_v0_over_cs - scaled ambient rupture velocity
%   rs_type        - 'slip' or 'aging' friction law
%   a_over_b       - friction parameter ratio
%   Delta_f0_over_b- scaled overstress
%   v0_over_cs     - scaled ambient sliding velocity
%   V_min          - approximate minimum of 𝒱 function
%
% Output:
%   vr_over_cs     - initial rupture velocity (v_r / c_s)

    % Define residual function handle using modular EoM_objective_func
    res_func = @(vr) EoM_objective_func(l_over_lb, vr, Param);

    % Bounds for rupture velocity
    vmin_over_cs = Param.V_min * Param.bar_v0_over_cs;       % small positive number
    vmax_over_cs = 1 - 1e-12;   % slightly below 1

    % Check if root exists in interval
    if res_func(vmin_over_cs) * res_func(vmax_over_cs) > 0
        error('No root exists in the given interval. Try a smaller initial crack length.');
    end

    % figure
    % hold on
    % x = linspace(vmin_over_cs,vmax_over_cs,1000);
    % for i = 1:length(x)
    %     resans(i) = res_func(x(i));
    % end
    % plot(x,resans,'o')

    % Solve for root using fzero
    options = optimset('Display','off');  % suppress output
    vr_over_cs = fzero(res_func, [vmin_over_cs, vmax_over_cs], options);

    % Optional: print
    fprintf('Initialized rupture velocity: v_r/c_s = %.6e\n', vr_over_cs);
end
