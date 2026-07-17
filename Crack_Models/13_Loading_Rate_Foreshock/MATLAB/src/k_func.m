function k = k_func(vr_over_cs)
    k = sqrt(1 - vr_over_cs.^2) ./ ellipe_safe(vr_over_cs);
end
