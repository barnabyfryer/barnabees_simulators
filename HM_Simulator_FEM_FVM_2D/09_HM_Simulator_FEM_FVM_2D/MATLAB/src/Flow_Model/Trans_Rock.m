function [Trans,dTransdP_LB,dTransdP_RT] = Trans_Rock(Flow,Gen,State)
%% - Rock Only Transmissibility
%Solves for the transimissibilities (non-fluid) in a 1D model. Constant
%grid block sizes. 
%Barnaby Fryer - EPFL - 29.03.17

%% - Find permeability
%Calculate perm and derivative
[State,~,dkx_dP,dky_dP,~] = Perm(Flow,State);                   %[N,1]

%Permeability [m^2]
kx = reshape(State.kx,Gen.Nx,Gen.Ny);                           %[N,1]
ky = reshape(State.ky,Gen.Nx,Gen.Ny);                           %[N,1]
%Derivative of permeability wrt pressure [m^2/Pa]
dkxdp = reshape(dkx_dP,Gen.Nx,Gen.Ny);                          %[N,1]
dkydp = reshape(dky_dP,Gen.Nx,Gen.Ny);                          %[N,1]

%% - Basic preparation
%Length of cell in x-direction [m]
dx = Gen.dx;                                                    %[1,1]
%Length of cell in y-direction [m]
dy = Gen.dy;                                                    %[1,1]       
%Length of cell in z-direction [m]
dz = Gen.Lz/1;                                                  %[1,1] 

%Cross-sectional area in x-direction [m^2]
Ax = dy*dz;                                                     %[1,1] 
%Cross-sectional area in y-direction [m^2]
Ay = dx*dz;                                                     %[1,1] 

%Initialize x-direction trans (of left of cell edge)
Trans.x = zeros(Gen.Nx+1,Gen.Ny);                               %[Nx+1,Ny]
dTransdP_LB.x = zeros(Gen.Nx+1,Gen.Ny);                         %[Nx+1,Ny]
dTransdP_RT.x = zeros(Gen.Nx+1,Gen.Ny);                         %[Nx+1,Ny]
%Initialize y-direction trans (of bottom of cell edge)
Trans.y = zeros(Gen.Nx,Gen.Ny+1);                               %[Nx,Ny+1]
dTransdP_LB.y = zeros(Gen.Nx,Gen.Ny+1);                         %[Nx,Ny+1]
dTransdP_RT.y = zeros(Gen.Nx,Gen.Ny+1);                         %[Nx,Ny+1]

%% - Do transmissibilities
%Find rock transmissibilities
Trans.x(2:Gen.Nx,:) = 2*kx(1:end-1,:).*kx(2:end,:)*Ax./...
    ((kx(1:end-1,:)+kx(2:end,:))*dx);                       %[Nx+1,Ny]
Trans.y(:,2:Gen.Ny) = 2*ky(:,1:end-1).*ky(:,2:end)*Ay./...
    ((ky(:,1:end-1)+ky(:,2:end))*dy);                       %[Nx,Ny+1]

%Find left edge
Trans.x(1,:) = kx(1,:)*Ax./(dx/2);                         	%[Nx+1,Ny]
%Find right edge
Trans.x(end,:) = kx(end,:)*Ax./(dx/2);                   	%[Nx+1,Ny]
%Find top edge
Trans.y(:,1) = ky(:,1)*Ay./(dy/2);                         	%[Nx,Ny+1]
%Find right edge
Trans.y(:,end) = ky(:,end)*Ay./(dy/2);                   	%[Nx,Ny+1]

%% - Do derivatives of transmissibilities
%Derivative of transmissibility wrt left or bottom cell
dTransdP_LB.x(2:Gen.Nx,:) = (2*Ax/dx) .* kx(2:end,:).^2 ./ (kx(1:end-1,:) + kx(2:end,:)).^2 .* dkxdp(1:end-1,:);%[Nx+1,Ny]
%Right most interface dependence on cell Nx
dTransdP_LB.x(Gen.Nx+1,:) = (2*Ax/dx) * dkxdp(end,:);

%Derivative of transmissibility wrt left or bottom cell
dTransdP_LB.y(:,2:Gen.Ny) = (2*Ay/dy) .* ky(:,2:end).^2 ./ (ky(:,1:end-1) + ky(:,2:end)).^2 .* dkydp(:,1:end-1);%[Nx,Ny+1]
%Right most interface dependence on cell Ny
dTransdP_LB.y(:,Gen.Ny+1) = (2*Ay/dy) * dkydp(:,end);

%Derivative of transmissibility wrt right or top cell
dTransdP_RT.x(2:Gen.Nx,:) = (2*Ax/dx) .* kx(1:end-1,:).^2 ./ (kx(1:end-1,:) + kx(2:end,:)).^2 .* dkxdp(2:end,:);%[Nx+1,Ny]
%Left most interface dependence on cell 1
dTransdP_RT.x(1,:) = (2*Ax/dx) * dkxdp(1,:);

%Derivative of transmissibility wrt right or top cell
dTransdP_RT.y(:,2:Gen.Ny) = (2*Ay/dy) .* ky(:,1:end-1).^2 ./ (ky(:,1:end-1) + ky(:,2:end)).^2 .* dkydp(:,2:end);%[Nx,Ny+1]
%Left most interface dependence on cell 1
dTransdP_RT.y(:,1) = (2*Ay/dy) * dkydp(:,1);

end

