function [Trans,dTransdP_LB,dTransdP_RT] = Trans_Rock(Flow,Gen,State)
%% - Rock Only Transmissibility
%Solves for the transimissibilities (non-fluid) in a 1D model. Constant
%grid block sizes. 
%Barnaby Fryer - EPFL - 29.03.17

%% - Find permeability
%Calculate perm and derivative
[State,dkdP] = Perm(Flow,State);                	        %[N,1]

%Permeability [m^2]
kx = reshape(State.kx,Gen.Nx,1);                            %[Nx,1]
%Derivative of permeability wrt pressure [m^2/Pa]
dkxdp = reshape(dkdP,Gen.Nx,1);                             %[Nx,1]

%% - Basic preparation
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
dTransdP_LB.x = zeros(Gen.Nx,1);                            %[Nx,1]
dTransdP_RT.x = zeros(Gen.Nx,1);                            %[Nx,1]

%% - Do transmissibilities
%Find rock transmissibilities
Trans.x(2:Gen.Nx,:) = 2*kx(1:end-1,:).*kx(2:end,:)*Ax./...
    ((kx(1:end-1,:)+kx(2:end,:))*dx);                       %[Nx+1,1]

%Find left edge
Trans.x(1,:) = kx(1,:)*Ax./(dx/2);                         	%[Nx+1,1]
%Find right edge
Trans.x(end,:) = kx(end,:)*Ax./(dx/2);                   	%[Nx+1,1]

%% - Do derivatives of transmissibilities
%Derivative of transmissibility wrt left or bottom cell
dTransdP_LB.x(2:Gen.Nx,:) = (2*Ax/dx) .* kx(2:end,:).^2 ./ (kx(1:end-1,:) + kx(2:end,:)).^2 .* dkxdp(1:end-1,:);%[Nx+1,1]
%Right most interface dependence on cell Nx
dTransdP_RT.x(Gen.Nx+1,:) = (2*Ax/dx) * dkxdp(end,:);

%Derivative of transmissibility wrt right or top cell
dTransdP_RT.x(2:Gen.Nx,:) = (2*Ax/dx) .* kx(1:end-1,:).^2 ./ (kx(1:end-1,:) + kx(2:end,:)).^2 .* dkxdp(2:end,:);%[Nx+1,1]
%Left most interface dependence on cell 1
dTransdP_RT.x(1,:) = (2*Ax/dx) * dkxdp(1,:);

end

