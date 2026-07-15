function zeta = Delta_fp_over_b(vr_over_cs, Param, kappa0)
    % DELTA_FP_OVER_B  Invert Δfp_over_b using G21's method

    %Find g
    g = g_func(vr_over_cs);

    % Value of 𝒱(ζ)
    V = vr_over_cs ./ Param.bar_v0_over_cs ./ g;  % use g directly

    % Initialize ζ_0
    zeta = log(V ./ kappa0);

    % Loop for iterative update
    n_iter = 2;
    for i = 1:n_iter
        zeta = log(V ./ chi_func(zeta, Param.rs_type, kappa0));
    end
end