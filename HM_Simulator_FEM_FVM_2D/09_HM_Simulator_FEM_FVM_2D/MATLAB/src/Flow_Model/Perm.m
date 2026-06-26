function [State,keff,dkx_dP,dky_dP,dkeff_dP] = Perm(Flow,State)

%Find each direction's permeability
State.kx = Flow.kx0 .* exp(Flow.ck*(State.P-Flow.kP0));
State.ky = Flow.ky0 .* exp(Flow.ck*(State.P-Flow.kP0));
%Add contribution of volumetric strain
State.kx = State.kx.*exp(-Flow.ckv * State.e_vol);
State.ky = State.ky.*exp(-Flow.ckv * State.e_vol);

%Find an effective permeability of the cell for use in wells
keff = sqrt(State.kx.*State.ky);

%Derivative wrt pressure
dkx_dP = Flow.ck .* State.kx;
dky_dP = Flow.ck .* State.ky;

%Derivative wrt pressure for effective permeability
dkeff_dP = 0.5*(State.kx.*State.ky).^(-0.5) .* (State.kx.*dky_dP + dkx_dP.*State.ky);

end

