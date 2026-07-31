function [phi,dphidp] = PhiCalc(Flow,Gen,State)

%% - Porosity
if Gen.Mandel == 1
    %Porosity and its derivative wrt to pressure after Coussy 2007; assume
    %initial e_vol = 0
    phi = -Gen.biot.*(State.e_vol + State.e_vol*0) + State.P.*(Gen.biot - Flow.phi0)./Flow.ks + Flow.phi0;
    %% - Derivative wrt Pressure
    dphidp = (Gen.biot - Flow.phi0)./Flow.ks;

else
    %Add contribution due to volumetric strain
    phi0 = Flow.phi0 - (1 - Flow.phi0).*State.e_vol;
    %Find porosity change considering pore pressure
    phi = phi0 .* exp(Flow.cphi*(State.P-Flow.phiP0));

    %% - Derivative wrt Pressure
    dphidp = Flow.cphi .* phi;
end

end

