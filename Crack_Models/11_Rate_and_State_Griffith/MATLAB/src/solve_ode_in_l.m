function [t_over_ts, l_over_lb, vr_over_cs, V_eff, reason] = solve_ode_in_l(l_ini_over_lb, l_fin_over_lb, Param)
% SOLVE_ODE_IN_L
% Solve for the time evolution of rupture along a crack using second-order ODE
%
% Inputs:
%   l_ini_over_lb   - initial crack length (l/l_b)
%   l_fin_over_lb   - final crack length (l/l_b)
%   Delta_T         - scaled hypocentral force
%   V_min           - approximate minimum of 𝒱 function
%   bar_v0_over_cs  - scaled ambient rupture velocity
%   a_over_b        - friction parameter ratio
%   Delta_f0_over_b - scaled overstress
%   v0_over_cs      - scaled ambient sliding velocity
%   rs_type         - 'slip' or 'aging' friction law
%
% Outputs:
%   t_over_ts       - dimensionless time array
%   l_over_lb       - crack length array
%   vr_over_cs      - rupture velocity array
%   reason          - stopping reason (0: v<threshold, 1: reached l_fin)

%% Initial conditions
t_ini_over_ts = 0;

% Find initial rupture velocity using initialize.m
vr_ini_over_cs = initialize(l_ini_over_lb, Param);

% Set initial ODE state: y = [t, dt/dl]'
y0 = [t_ini_over_ts; 1 / vr_ini_over_cs];

%% Nested function: ODE RHS
    function dy = ode_rhs(l, y)
        %Guard rail to stop vr going past cs
        vr = 1 / y(2);
        if vr >= 1
            % Freeze dynamics to avoid invalid math, exits function
            dy = [y(2); 0];
            return
        end


        t_val = y(1);
        dt_dl = y(2);
        vr = 1 / dt_dl;

        % Compute derivatives using modular functions
        num = EoM_dl(l, vr, Param);
        den = EoM_dv(l, vr, Param);

        d2t_dl2 = (1 / vr^2) * (num / den);
        dy = [dt_dl; d2t_dl2];

        if rand < 1e-1
            fprintf('l = %.3e, v = %.3e\n', l, vr);
        end
    end

%% Event function to stop if rupture velocity drops below or above threshold
v_thres = Param.V_min * Param.bar_v0_over_cs;
%Combined function with both conditions
    function [value, isterminal, direction] = combined_events(l, y)
        [v1, t1, d1] = event_stop(l, y);
        [v2, t2, d2] = event_supersonic(l, y);
        value = [v1; v2];
        isterminal = [t1; t2];
        direction = [d1; d2];
    end
%Condition of vr < cs
    function [value, isterminal, direction] = event_supersonic(l, y)
        vr = 1 / y(2);
        value = 1 - vr;     % stop when vr >= 1
        isterminal = 1;
        direction = -1;
    end
%Condition of rupture halting
    function [value, isterminal, direction] = event_stop(l, y)
        value = 1 / y(2) - v_thres;  % stop when 1/y(2) < threshold
        isterminal = 1;
        direction = -1;  % negative crossing only
    end

%% Solve ODE using ode45
options = odeset('RelTol',Param.RelTol, 'AbsTol',Param.AbsTol,'Events', @combined_events);
% l_span = logspace(log10(l_ini_over_lb), log10(l_fin_over_lb), 20);

[l_sol, y_sol, te, ye, ie] = ode45(@ode_rhs, [l_ini_over_lb, l_fin_over_lb], y0, options);
% [l_sol, y_sol, te, ye, ie] = ode45(@ode_rhs, l_span, y0, options);

%% Extract outputs
l_over_lb = l_sol;
t_over_ts = y_sol(:,1);
vr_over_cs = 1 ./ y_sol(:,2);

%Solve for effective slip velocity
V_eff = zeros(size(l_over_lb));
for i = 1:length(l_over_lb)
    V_eff(i) = Veff_func(l_over_lb(i), vr_over_cs(i), Param);
end

%% Determine stopping reason
if ~isempty(ie)
    if ie(end) == 1
        reason = 0;  % rupture died
    elseif ie(end) == 2
        reason = 3;  % supersonic nucleation
    end
elseif l_over_lb(end) >= l_fin_over_lb
    reason = 1; % reached final crack length
else
    reason = 2; % integration incomplete
    warning('ODE solver did not reach final crack length.');
end

end




