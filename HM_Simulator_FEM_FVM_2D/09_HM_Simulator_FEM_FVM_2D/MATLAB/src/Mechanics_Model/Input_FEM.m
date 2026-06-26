function [Gen, Pos, State, Storage] = Input_FEM(Flow, Gen, State, Storage)


%% - General Inputs
%Young's Modulus [Pa]
Gen.E = 1e9;                                           %[1,1]
%Poisson's Ratio [-]
Gen.v = .3;                                             %[1,1]
%Biot's coefficient [-]
Gen.biot = 1;                                           %[1,1]
%Select plane stress (stress) or plane strain (strain)
Gen.plane = "strain";

%% - Basic Calculations
% Cell center coordinates
xc = linspace(Gen.dx/2, Gen.Lx-Gen.dx/2, Gen.Nx);
yc = linspace(Gen.dy/2, Gen.Ly-Gen.dy/2, Gen.Ny);

% Mesh of cell centers
[xc,yc] = meshgrid(xc,yc);                        
Storage.x = reshape(xc',1,Gen.Nx*Gen.Ny);               %[1,N]
Storage.y = reshape(yc',1,Gen.Nx*Gen.Ny);               %[1,N]

%Element edge lengths
Gen.dx = Gen.Lx/Gen.Nx;                                 %[1,1]
Gen.dy = Gen.Ly/Gen.Ny;                                 %[1,1]

%Number of total elements
Gen.Ne = Gen.Nx*Gen.Ny;                                 %[1,1]
%Total number of nodes
Gen.Nn = (Gen.Nx+1)*(Gen.Ny+1);                         %[1,1]
%Matrix with the nodes making up each element
Gen.Ref = zeros(Gen.Ne,Gen.Nn);                         %[Ne,Nn]
%Loop over all elements
for i = 1:Gen.Ne
    %Find the row that element is in
    row = ceil(i/Gen.Nx);                               %[1,1]
    %Assign bottom two nodes to be in that element
    Gen.Ref(i,i+row-1:i+row) = 1;                       %[Ne,Nn]
    %Assign the top two nodes to be in that element
    Gen.Ref(i,i+row+Gen.Nx:i+row+Gen.Nx+1) = 1;         %[Ne,Nn]
end

%Original node locations
Pos.x = ones(Gen.Nx+1,Gen.Ny+1).*linspace(0,Gen.Lx,Gen.Nx+1)';  %[Nx+1,Ny+1]
Pos.y = ones(Gen.Nx+1,Gen.Ny+1).*linspace(0,Gen.Ly,Gen.Ny+1);   %[Nx+1,Ny+1]

%Elements are counted as follows:
%Reservoir:             In matrix: 
%13 14 15 16            (1,4) (2,4) (3,4) (4,4) 
%09 10 11 12            (1,3) (2,3) (3,3) (4,3) 
%05 06 07 08            (1,2) (2,2) (3,2) (4,2) 
%01 02 03 04            (1,1) (2,1) (3,1) (4,1) 

%Nodes are counted as follows:
%Reservoir:             In matrix: 
%21 22 23 14 25           (1,5) (2,5) (3,5) (4,5) (5,5) 
%16 17 18 19 20           (1,4) (2,4) (3,4) (4,4) (5,4) 
%11 12 13 14 15           (1,3) (2,3) (3,3) (4,3) (5,3) 
%06 07 08 09 10           (1,2) (2,2) (3,2) (4,2) (5,2) 
%01 02 03 04 05           (1,1) (2,1) (3,1) (4,1) (5,1)

%% - Fixed displacement nodes
% %Fix x-direction
% %Fix left nodes
Gen.nodesx(1:Gen.Ny+1,1) = (1:Gen.Nx+1:Gen.Nn-Gen.Nx)';
% %Fix right nodes
% Gen.nodesx(Gen.Ny+2:2*Gen.Ny+2) = (Gen.Nx+1:Gen.Nx+1:Gen.Nn)';
%  
% %Fix y-directions
% %Fix bottom nodes
Gen.nodesy(1:Gen.Nx+1,1) = (1:Gen.Nx+1)';
% %Fix top nodes
% Gen.nodesy(Gen.Nx+2:2*Gen.Nx+2) = (Gen.Nn-Gen.Nx:Gen.Nn)';

%Fix middle node
% Gen.nodesx = floor((Gen.Ny+1)/2)*(Gen.Nx+1) - floor((Gen.Nx+1)/2);
% Gen.nodesx(2) = Gen.nodesx(1) + 1;
% Gen.nodesx(3) = Gen.nodesx(1) + Gen.Nx + 1;
% Gen.nodesx(4) = Gen.nodesx(2) + Gen.Nx + 1;
% Gen.nodesy = Gen.nodesx;

%% - Initialize phi and e_vol
%Initialize volumetric strain   
State.e_vol = zeros(Gen.Nx*Gen.Ny,1);                           %[N,1]
%Find the e_vol at initial conditions
[State] = M_Simulator_FEM_2D(Gen,Pos,State);

%Initialize permeability
[State,~,~,~,~] = Perm(Flow,State);                             %[N,1]
%Initialize porosity
[State.phi,~] = PhiCalc(Flow,State);                	        %[N,1]

%% - Storage matrices
%Number of points to store
TStore = min([200,floor(Gen.tf/Gen.tstep)]);                    %[1,1]
%Get storage times
Storage.TStorage = 0:Gen.tf/TStore:Gen.tf;                      %[1,TStore]
Storage.TStorage = floor(Storage.TStorage/Gen.tstep)*Gen.tstep; %[1,TStore]
%Predefine memory for saved pressures
Storage.P = zeros(TStore,Gen.Nx*Gen.Ny);                        %[TStore,N]
%Store initial pressure
Storage.P(1,:) = State.P;                                       %[TStore,N]

%Predefine memory for saved permeabilities
Storage.kx = zeros(TStore,Gen.Nx*Gen.Ny);                       %[TStore,N]
Storage.ky = zeros(TStore,Gen.Nx*Gen.Ny);                       %[TStore,N]
%Store initial pressure
Storage.kx(1,:) = State.kx;                                     %[TStore,N]
Storage.ky(1,:) = State.kx;                                     %[TStore,N]

%Predefine memory for saved porosities
Storage.phi = zeros(TStore,Gen.Nx*Gen.Ny);                      %[TStore,N]
%Store initial pressure
Storage.phi(1,:) = State.phi;                                   %[TStore,N]

%Predefine memory for saved flux
Storage.flux = zeros(TStore,Gen.Nx*Gen.Ny);                     %[TStore,N]

%Predefine memory for saved volumetric strain
Storage.e_vol = zeros(TStore,Gen.Nx*Gen.Ny);                    %[TStore,N]
%Store initial volumetric strain
Storage.e_vol(1,:) = State.e_vol;                               %[TStore,N]

%Predefine memory for saved stresses
Storage.Sig_xx = zeros(TStore,Gen.Nx*Gen.Ny);                   %[TStore,N]
Storage.Sig_yy = zeros(TStore,Gen.Nx*Gen.Ny);                   %[TStore,N]
Storage.Sig_xy = zeros(TStore,Gen.Nx*Gen.Ny);                   %[TStore,N]
%Store initial volumetric strain
Storage.Sig_xx(1,:) = State.Sig_xx;                             %[TStore,N]
Storage.Sig_yy(1,:) = State.Sig_xx;                             %[TStore,N]
Storage.Sig_xy(1,:) = State.Sig_xx;                             %[TStore,N]

%% - Plotting inputs
Plotting.lwidth_1col = 0.75;
Plotting.Position_1col_matrix = [2.2 1.8 6 4.5];
Plotting.fsize_1col = 7;

end

