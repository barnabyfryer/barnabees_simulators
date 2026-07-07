function E = ellipe_safe(v)
    m = min(max(1 - v.^2, 0.0), 1.0 - 1e-14);
    [~, E] = ellipke(m);  % second output is complete elliptic integral of 2nd kind
end