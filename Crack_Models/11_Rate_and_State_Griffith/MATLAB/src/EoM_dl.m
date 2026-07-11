function dFdl = EoM_dl(l, v, Param)

    h = 1e-6 * max(abs(l), realmin);

    dFdl = (EoM_objective_func(l+h, v, Param) ...
          - EoM_objective_func(l-h, v, Param)) / (2*h);
end
