function [Trans] = Trans_Rock(Flow,Gen)
%% - Rock Only Transmissibility
%Solves for the transimissibilities (non-fluid) in a 1D model. Constant
%grid block sizes. Homogeneous.
%Barnaby Fryer - EPFL - 29.03.17

%Permeability [m^2]
kx = reshape(Flow.kx,Gen.Nx,1);                             %[Nx,1]

%Length of cell in x-direction [m]
dx = Gen.Lx/Gen.Nx;                                         %[1,1]
%Length of cell in y-direction [m]
dy = Gen.Ly/1;                                              %[1,1]       
%Length of cell in z-direction [m]
dz = Gen.Lz/1;                                              %[1,1] 

%Cross-sectional area in x-direction [m^2]
Ax = dy*dz;                                                 %[1,1] 

%Initialize x-direction trans (of left of cell edge)
Trans.x = zeros(Gen.Nx+1,1);                                %[Nx+1,1]

%Find rock transmissibilities
Trans.x(2:Gen.Nx,:) = 2*kx(1:end-1,:).*kx(2:end,:)*Ax./...
    ((kx(1:end-1,:)+kx(2:end,:))*dx);                       %[Nx+1,1]

%Find left edge
Trans.x(1,:) = kx(1,:)*Ax./(dx/2);                         	%[Nx+1,1]
%Find right edge
Trans.x(end,:) = kx(end,:)*Ax./(dx/2);                   	%[Nx+1,1]

end

