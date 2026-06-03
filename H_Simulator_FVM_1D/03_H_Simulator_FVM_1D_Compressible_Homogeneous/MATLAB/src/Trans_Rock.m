function [Trans] = Trans_Rock(Flow,Gen)
%% - Rock Only Transmissibility
%Solves for the transimissibilities (non-fluid) in a 1D model. Constant
%grid block sizes. Homogeneous.
%Barnaby Fryer - EPFL - 29.03.17

%Permeability [m^2]
kx = Flow.kx;                                               %[Nx,1]

%Length of cell in x-direction [m]
dx = Gen.Lx/Gen.Nx;                                         %[1,1]
%Length of cell in y-direction [m]
dy = Gen.Ly/1;                                              %[1,1]       
%Length of cell in z-direction [m]
dz = Gen.Lz/1;                                              %[1,1] 

%Cross-sectional area in x-direction [m^2]
Ax = dy*dz;                                                 %[1,1] 

%Find rock transmissibilities
Trans.x = kx*Ax./dx;                                        %[1,1]

end

