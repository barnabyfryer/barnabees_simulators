function [t,l,v,df0_over_b] = implicit_continuation(...
        t0,...
        l0,...
        v0,...
        Param)

%Switch to an implicit solver in the case that the manifold becomes too
%stiff. Now we solve for l as a function of v (instead of v as function of
%l as before). 

%Number of discretization points in v
N = 100;

%Go up to what velocity doing this
vmax = 0.95;

%Build velocities to interrogate
v = logspace(log10(v0),log10(vmax),N);

%Predefine memory for t and l
l = zeros(size(v));
t = zeros(size(v));
df0_over_b = zeros(size(v));

%Set initial condition where explicit solver failed
l(1)=l0;
t(1)=t0;
df0_over_b(1) = Param.Delta_f0_over_b_in + Param.Loading_rate * t0 * Param.v0_over_cs;

for i=2:N
    t_new = t(i-1);
    error = 1;
    while error > Param.RelTol
        %Update overstress
        Param.Delta_f0_over_b = Param.Delta_f0_over_b_in + Param.Loading_rate * t_new * Param.v0_over_cs;
        %Set up function solver for l
        fun = @(L) EoM_objective_func(L,v(i),Param);
        %Solve for new l using previous l
        l_old = fzero(fun,l(i-1));
        
        %Find average velocity of this step
        vavg = 0.5*(v(i) + v(i-1));
        %Time of step
        dt = (l_old - l(i-1))/vavg;
        %Update time
        t_new = t(i-1) + dt;
        %Update overstress
        Param.Delta_f0_over_b = Param.Delta_f0_over_b_in + Param.Loading_rate * t_new * Param.v0_over_cs;
        %Remake function solver for l
        fun = @(L) EoM_objective_func(L,v(i),Param);
        %New l
        l_new = fzero(fun,l_old);
        %Error 
        error = abs(l_old - l_new)/max(abs(l_new),1);
    end
    %Save time
    t(i) = t_new;
    %Save crack length
    l(i) = l_new;
    %Save overstress
    df0_over_b(i) = Param.Delta_f0_over_b;
end

end