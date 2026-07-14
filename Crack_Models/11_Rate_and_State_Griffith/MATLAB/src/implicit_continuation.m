function [t,l,v] = implicit_continuation(...
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

%Set initial condition where explicit solver failed
l(1)=l0;
t(1)=t0;

for i=2:N
    %Set up function solver for l
    fun = @(L) EoM_objective_func(L,v(i),Param);
    %Solve for new l using previous l
    l(i)=fzero(fun,l(i-1));

end

vavg = 0.5*(v(1:end-1)+v(2:end));

dt = diff(l)./vavg;

t(2:end)=t0+cumsum(dt);

end