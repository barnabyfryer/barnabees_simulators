function dFdv = EoM_dv(l, v, Param)

    h = 1e-6 * max(abs(v), realmin);

    if v-h <= 0
        dFdv = (EoM_objective_func(l,v+h,Param) ...
              - EoM_objective_func(l,v,Param))/h;
    else
        dFdv = (EoM_objective_func(l,v+h,Param) ...
              - EoM_objective_func(l,v-h,Param))/(2*h);
    end

end
