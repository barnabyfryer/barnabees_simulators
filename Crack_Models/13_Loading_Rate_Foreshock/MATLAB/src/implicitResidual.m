function [res] = implicitResidual(l_over_lb_guess, l_over_lb_old, t_old, vavg, vr_over_cs_new, Param)
%Update time
t = t_old + (l_over_lb_guess-l_over_lb_old)/vavg;
%Update overstress
Param.Delta_f0_over_b = ...
    Param.Delta_f0_over_b_in + Param.Loading_rate*(t+t_old)/2*Param.v0_over_cs;
%Find residual
res = EoM_objective_func(l_over_lb_guess,vr_over_cs_new,Param);
end