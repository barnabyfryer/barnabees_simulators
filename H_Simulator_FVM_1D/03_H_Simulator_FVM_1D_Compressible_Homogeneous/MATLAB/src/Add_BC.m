function [BC,dBC] = Add_BC(dRhodP,Gen,Rho,State,Trans)

%Predefine memory
BC = zeros(Gen.Nx,1);                                           %[N,1]
dBC = zeros(Gen.Nx,1);                                          %[N,1]

%Find "source" terms from boundary conditions
BC(1) = (Trans.x/2) * Rho(1) * (State.P(1)  - Gen.PL);          %[N,1]
BC(end) = (Trans.x/2) * Rho(end) * (State.P(end)- Gen.PR);      %[N,1]


dBC(1) = (Trans.x/2) * (Rho(1) + dRhodP(1)*(State.P(1) - Gen.PL));
dBC(end) = (Trans.x/2) * (Rho(end) + dRhodP(end)*(State.P(end) - Gen.PR));

end