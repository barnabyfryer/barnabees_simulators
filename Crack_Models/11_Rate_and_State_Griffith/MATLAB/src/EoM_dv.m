function dFdv = EoM_dv(l, v, Param)

    h = 1e-6 * max(abs(v), realmin);

    dFdv = (EoM_objective_func(l, v+h, Param) ...
          - EoM_objective_func(l, v-h, Param)) / (2*h);
end
