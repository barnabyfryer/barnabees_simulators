function Kc = Kc_func(vr_over_cs, Param)
    Gc = Gc_func(vr_over_cs, Param);
    g = g_func(vr_over_cs);
    k = k_func(vr_over_cs);
    Kc = sqrt(2 .* Gc .* g ./ k.^2);
end
