function [State,dkdP] = Perm(Flow,State)

State.kx = Flow.kx0 .* exp(Flow.ck*(State.P-Flow.kP0));

dkdP = Flow.ck .* State.kx;

end

