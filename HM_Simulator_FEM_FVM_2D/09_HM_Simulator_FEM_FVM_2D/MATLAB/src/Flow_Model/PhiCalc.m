function [phi,dphidp] = PhiCalc(Flow,State)

%% - Porosity
%Add contribution due to volumetric strain
phi0 = Flow.phi0 - (1 - Flow.phi0).*State.e_vol;
%Find porosity change considering pore pressure
phi = phi0 .* exp(Flow.cphi*(State.P-Flow.phiP0));

%% - Derivative wrt Pressure
dphidp = Flow.cphi .* phi;

end

