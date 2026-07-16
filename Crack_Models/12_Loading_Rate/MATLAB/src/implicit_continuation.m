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
    %Average velocity
    vavg = 0.5*(v(i) + v(i-1));
    %Implicit residual which computes with updated df0
    fun = @(L) implicitResidual(L, ...
        l(i-1), ...
        t(i-1), ...
        vavg, ...
        v(i), ...
        Param);

    %Evaluate crack length
    l(i) = fzero(fun,l(i-1));
    %Update time
    t(i) = t(i-1) + (l(i)-l(i-1))/vavg;
    %Save overstress
    df0_over_b(i) = Param.Delta_f0_over_b_in + Param.Loading_rate * t(i) * Param.v0_over_cs;
end

end