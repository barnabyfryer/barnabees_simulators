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
[vr_ini_over_cs, l_ini_over_lb] = initialize(l_ini_over_lb, Param);

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


        t_val = y(1); %t/ts * cs/v0
        dt_dl = y(2);
        vr = 1 / dt_dl;

        % Compute derivatives using modular functions
        num = EoM_dl(l, vr, Param);
        den = EoM_dv(l, vr, Param);

        %dt_dl = 1/v
        %d2t_dl2 = (1/v^2)(dF/dl)/(dF/dv)

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
% options = odeset('RelTol',Param.RelTol, 'AbsTol',Param.AbsTol,'Events', @combined_events,'InitialStep',1e-6,'MaxStep',1e-2);
options = odeset('RelTol',Param.RelTol, 'AbsTol',Param.AbsTol,'Events', @combined_events,'InitialStep',1e-6);

[l_sol, y_sol, te, ye, ie] = ode45(@ode_rhs, [l_ini_over_lb, l_fin_over_lb], y0, options);

%% Extract outputs
l_over_lb = l_sol;
t_over_ts = y_sol(:,1);
vr_over_cs = 1 ./ y_sol(:,2);

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
    warning('ODE explicit solver did not reach final crack length.');

        [t_imp,l_imp,v_imp] = implicit_continuation( ...
            t_over_ts(end),...
            l_over_lb(end),...
            vr_over_cs(end),...
            Param);

        %Get solved indices
		idx = find(isfinite(l_imp), 1, 'last');
		%Get solved values 
		l0 = l_imp(idx);
		t0 = t_imp(idx);
		v0 = v_imp(idx);

		%Get restart conditions
		y0_restart = [t0; 1/v0];

		%Start using explicit method to finish
		[l_sol, y_sol, te, ye, ie] = ode45(@ode_rhs, [l0, l_fin_over_lb], y0_restart, options);

		%Extract end values for second explicit section
		l_over_lb_end = l_sol;
		t_over_ts_end = y_sol(:,1);
		vr_over_cs_end = 1 ./ y_sol(:,2);

		t_over_ts = [t_over_ts; t_imp(1:idx)'; t_over_ts_end];

		l_over_lb = [l_over_lb; l_imp(1:idx)'; l_over_lb_end];

		vr_over_cs = [vr_over_cs; v_imp(1:idx)'; vr_over_cs_end];
end

%% - Solve for effective slip velocity
V_eff = zeros(size(l_over_lb));
for i = 1:length(l_over_lb)
    V_eff(i) = Veff_func(l_over_lb(i), vr_over_cs(i), Param);
end

end




