function [Rho,dRhodP] = Density(Flow,P)
%% - Density Calculator
%Finds the density and the derivative of density with pressure using the
%definition (EOS) of a slightly compressible fluid. The equations is
%vectorized so it takes a vector (for pressure) and spits out a vector.
%This is for slightly compressible fluids
%Barnaby Fryer - EPFL - 29.03.17


%Compressibility [1/Pa]
cf = Flow.cf;                                           %[1,1]
%Reference density [kg/m^3]
Rho0 = Flow.Rho0;                                       %[1,1]
%Reference pressure [Pa]
RhoP = Flow.RhoP;                                       %[1,1]

%Compute density
Rho = Rho0 * exp(cf*(P-RhoP));                          %[N,1]
%Compute derivative wrt pressure
dRhodP = Rho*cf;                                        %[N,1]
end

