%% - Input function
function [Flow, Gen, Plotting, State, Storage, Wells] = InputFlowData()
%% - General Parameters
%Final time [sec]
Gen.tf = 200;                                                   %[1,1]    
%Time step [sec]
Gen.tstep = 5;                                                  %[1,1]    
%Tolerance [-]
Gen.tol = 1e-5;                                                 %[1,1]   
%Tolerance of sequential coupling [-]
Gen.tol_all = 1e-3;                                             %[1,1]
%Number of cells in x-direction [-]
Gen.Nx = 41;                                                    %[1,1]  
%Number of cells in y-direction [-]
Gen.Ny = 41;                                                    %[1,1]  
%Reservoir length in x-direction [m]
Gen.Lx = 10;                                                    %[1,1]   
%Reservoir length in y-direction [m]
Gen.Ly = 10;                                                    %[1,1]    
%Reservoir height [m]
Gen.Lz = 1;                                                     %[1,1]  
%Initial pressure in reservoir [Pa]
Gen.Pi = 1e7;

%% - Basic Calculations
%Total number of elements
Gen.N = Gen.Nx*Gen.Ny;                                          %[1,1]
%Element edge lengths
Gen.dx = Gen.Lx/Gen.Nx;                                         %[1,1]
%Element edge lengths
Gen.dy = Gen.Ly/Gen.Ny;                                         %[1,1]
%Element edge lengths
Gen.dz = 1;                                                     %[1,1]

% Cell center coordinates
xc = linspace(Gen.dx/2, Gen.Lx-Gen.dx/2, Gen.Nx);
yc = linspace(Gen.dy/2, Gen.Ly-Gen.dy/2, Gen.Ny);

% Mesh of cell centers
[xc,yc] = meshgrid(xc,yc);                        
Storage.x = reshape(xc',1,Gen.Nx*Gen.Ny);                       %[1,N]
Storage.y = reshape(yc',1,Gen.Nx*Gen.Ny);                       %[1,N]

%% - Flow Model
%x-direction
%Permeability, of left edges [m^2]
kx_left = 1e-12;                                                %[1,1] 
%Permeability, of right zone [m^2]
kx_right = 1e-12;                                               %[1,1] 
%Permeability left edge bands length [m]
Lx_k = 2;                                                       %[1,1] 
%Permeability [m^2] 
Flow.kx0 = kx_right*ones(Gen.Nx*Gen.Ny,1);                      %[1,1]
Flow.kx0(Storage.x < Lx_k) = kx_left;

%y-direction
%Permeability, of left edges [m^2]
ky_left = 1e-12;                                                %[1,1] 
%Permeability, of right zone [m^2]
ky_right = 1e-12;                                               %[1,1] 
%Permeability left edge bands length [m]
Ly_k = 2;                                                       %[1,1] 
%Permeability [m^2] 
Flow.ky0 = ky_right*ones(Gen.Nx*Gen.Ny,1);                      %[1,1]
Flow.ky0(Storage.y < Ly_k) = ky_left;

%"Compressibility" of permeability
Flow.ck = 1e-8;                                                 %[1,1] 
%"Compressibility" of due to volumetric strain
Flow.ckv = 1e-8;                                                %[1,1] 
%Reference pressure [Pa]
Flow.kP0 = 1e5;                                                 %[1,1]  

%Porosity, of left edges [-]
phi_left = 0.2;                                                 %[1,1] 
%Porosity, of right zone [-]
phi_right = 0.2;                                                %[1,1] 
%Porosity left edge bands length [m]
L_phi = 2;                                                      %[1,1] 
%Porosity [-]
Flow.phi0 = phi_right*ones(Gen.Nx*Gen.Ny,1);                    %[1,1]  
Flow.phi0(Storage.x < L_phi) = phi_left;

%"Compressibility" of porosity
Flow.cphi = 1e-9;                                               %[1,1] 
%Reference pressure [Pa]
Flow.phiP0 = 1e5;                                               %[1,1] 


%Fluid compressibility [1/Pa]
Flow.cf = 1e-8;                                                 %[1,1]          
%Fluid viscosity [Pa sec]
Flow.muf = .1;                                                  %[1,1]           
%Reference density [kg/m^3]
Flow.Rho0 = 1000;                                               %[1,1]         
%Reference pressure [Pa]
Flow.RhoP = 1e5;                                                %[1,1]  

%% - Wells
% Constant pressure wells
%Row vector of well pressures [Pa]
Wells.P = [3.5e7];                                              %[1,Nwells]
%Well indexes [m]
Wells.WI = [1000];                                              %[1,Nwells]
% Well locations [m]
Wells.xP = [5];                                                 %[1,Nwells]
Wells.yP = [5];                                                 %[1,Nwells]

    
%Find the cells of these wells
Wells.Loc_P = zeros(size(Wells.xP));
for i = 1:length(Wells.xP)
    [~, Wells.Loc_P(i)] = min((Storage.x - Wells.xP(i)).^2 + (Storage.y - Wells.yP(i)).^2);
end

%Constant rate wells (always keep at least a zero contribution in one cell
%Define constant rate [kg/sec]
Wells.Q = 0;                                                    %[1,Nwells_Q]
Wells.xQ = 5;                                                   %[1,Nwells_Q]
Wells.yQ = 5;                                                   %[1,Nwells_Q]

%Find the cells of these wells
Wells.Loc_Q = zeros(size(Wells.xQ));
for i = 1:length(Wells.xQ)
    [~, Wells.Loc_Q(i)] = min((Storage.x - Wells.xQ(i)).^2 + (Storage.y - Wells.yQ(i)).^2);
end


%% - Initialize state
%Initialize time
State.t = 0;                                                    %[1,1]
%Initialize pressure
State.P = zeros(Gen.Nx*Gen.Ny,1) + Gen.Pi;                      %[N,1]

%Time step counter
State.step = 1;                                                 %[1,1]

%% - Plotting inputs
Plotting.lwidth_1col = 0.75;
Plotting.Position_1col_matrix = [2.2 1.8 6 4.5];
Plotting.fsize_1col = 7;




end