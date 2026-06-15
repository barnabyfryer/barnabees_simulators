function [phi,dphidp] = PhiCalc(Flow,State)

%% - Porosity
phi = Flow.phi0 .* exp(Flow.cphi*(State.P-Flow.phiP0));

%% - Derivative wrt Pressure
dphidp = Flow.cphi .* phi;

end

