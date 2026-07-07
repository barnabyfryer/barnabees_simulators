function dFdl = EoM_dl(l, v, Param)
    h = 1e-8;
    dFdl = (EoM_objective_func(l+h, v, Param) ...
          - EoM_objective_func(l-h, v, Param)) / (2*h);
end
