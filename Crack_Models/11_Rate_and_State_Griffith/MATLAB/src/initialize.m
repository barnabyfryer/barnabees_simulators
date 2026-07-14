
function [vr_over_cs, l_over_lb] = initialize(l_over_lb, Param)
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


% Solve for root using fzero
options = optimset('Display','off');  % suppress output
vr_over_cs = fzero(res_func, [vmin_over_cs, vmax_over_cs], options);

%Find residual of initialization
F = abs(EoM_objective_func(l_over_lb,vr_over_cs,Param));
vnew = vr_over_cs;

%Check to see if initial crack length is too small
max_it = 100;
it = 0;
%Updates to initial crack length if too small
l_in_updates_allowed = 3;
l_updates_count = 0;
while F > 1e-15
    %Use new velocity
    v = vnew;
    %Find residual
    F = EoM_objective_func(l_over_lb,v,Param);
    %Find residual dependence on vr
    Fv = EoM_dv(l_over_lb, v, Param);
    %Update velocity
    vnew = v - F/Fv;
    %Re-evaluate residual
    F = abs(EoM_objective_func(l_over_lb,vnew,Param));
    %Update number of iterations
    it = it + 1;
    if it >= max_it
        if l_updates_count < l_in_updates_allowed
            %Count updates to initial crack length
            l_updates_count = l_updates_count + 1;
            %Reset iteration count
            it = 0;
            %Increase initial crack length to put solution on correct
            %manifold
            l_over_lb = l_over_lb*10;
            %Let the user know we're changing their crack length
            fprintf('Increasing initial crack length to: l/l_b = %.6e\n', l_over_lb);
            %Find initial rupture velocity for this crack length
            vnew = fzero(res_func, [vmin_over_cs, vmax_over_cs], options);
        else
            %Kick error, crack length likely too small
            error(['Initialization failed after %d Newton iterations. ' ...
                'Residual = %.3e. Initial crack length l/l_b = %.3e is likely too small.'], ...
                it, F, l_over_lb);
        end
    end

end

%Store final value of vr
vr_over_cs = vnew;

%Print
fprintf('Initialized rupture velocity: v_r/c_s = %.6e\n', vr_over_cs);
end
