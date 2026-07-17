function Gc = Gc_func(vr_over_cs, Param)
% GC_FUNC  Scaled fracture energy from Garagash (2021), Eq. (2.12)

    % Prefactor κ0
    switch Param.rs_type
        case "slip"
            kappa0 = 0.838;
        case "aging"
            kappa0 = 0.916;
        otherwise
            error("Unknown rate-and-state law.")
    end

    % Compute ζ = Δfp / b
    zeta = Delta_fp_over_b(vr_over_cs, Param, kappa0);

    % Apply constitutive law
    switch Param.rs_type
        case "slip"
            Gc = zeta;

        case "aging"
            Gc = -dilog_func(1 - exp(zeta));
    end
end
