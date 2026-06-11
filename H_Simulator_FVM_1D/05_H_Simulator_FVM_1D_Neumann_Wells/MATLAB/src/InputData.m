%% - Input function
function [Flow, Gen, Plotting, State, Storage, Wells] = InputData()
%% - General Parameters
%Final time [sec]
Gen.tf = 10000;                                                 %[1,1]    
%Time step [sec]
Gen.tstep = 10;                                                %[1,1]    
%Tolerance [-]
Gen.tol = 1e-5;                                                 %[1,1]    
%Number of cells in x-direction [-]
Gen.Nx = 500;                                                  %[1,1]  
%Reservoir length in x-direction [m]
Gen.Lx = 10;                                                    %[1,1]   
%Reservoir length in y-direction [m]
Gen.Ly = 1;                                                     %[1,1]    
%Reservoir height [m]
Gen.Lz = 1;                                                     %[1,1]  

%% - Basic Calculations
%Element edge lengths
Gen.dx = Gen.Lx/Gen.Nx;                                         %[1,1]
%Cell locations
Storage.x = linspace(Gen.dx/2,Gen.Lx-Gen.dx/2,Gen.Nx);          %[N,1]

%% - Flow Model
%Permeability, of left edges [m^2]
k_left = 1e-12;                                                  %[1,1] 
%Permeability, of right zone [m^2]
k_right = 1e-13;                                                   %[1,1] 
%Permeability left edge bands length [m]
L_k = 2;                                                        %[1,1] 
%Permeability [m^2] 
Flow.kx = k_right*ones(Gen.Nx,1);                                  %[1,1]
Flow.kx(Storage.x < L_k) = k_left;
%Porosity, of left edges [-]
phi_left = 0.3;                                                  %[1,1] 
%Porosity, of right zone [-]
phi_right = 0.2;                                                   %[1,1] 
%Porosity left edge bands length [m]
L_phi = 2;                                                      %[1,1] 
%Porosity [-]
Flow.phi = phi_right*ones(Gen.Nx,1);                               %[1,1]  
Flow.phi(Storage.x < L_phi) = phi_left;
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
Wells.P = [5e7 1e7];                                            %[1,Nwells]
%Well indexes [m^3]
Wells.WI = [1 1];                                               %[1,Nwells]
% Well locations [m]
Wells.xP = [10 0];                                              %[1,Nwells]
    
%Find the cells of these wells
Wells.Loc_P = zeros(size(Wells.xP));
for i = 1:length(Wells.xP)
    [~, Wells.Loc_P(i)] = min(abs(Storage.x - Wells.xP(i)));
end

%Constant rate wells (always keep at least a zero contribution in one cell
%Define constant rate [kg/sec]
Wells.Q = 1;                                                    %[1,Nwells_Q]
Wells.xQ = 4;                                                   %[1,Nwells_Q]

%Find the cells of these wells
Wells.Loc_Q = zeros(size(Wells.xQ));
for i = 1:length(Wells.xQ)
    [~, Wells.Loc_Q(i)] = min(abs(Storage.x - Wells.xQ(i)));
end


%% - Initialize state
%Initialize time
State.t = 0;                                                    %[1,1]
%Initialize pressure
State.P = zeros(Gen.Nx,1) + 1e5;                                %[N,1]
%Time step counter
State.step = 1;                                                 %[1,1]

%% - Storage matrices
%Number of points to store
TStore = 5;                                                     %[1,1]
%Get storage times
Storage.TStorage = 0:Gen.tf/TStore:Gen.tf;                      %[1,TStore]
Storage.TStorage = floor(Storage.TStorage/Gen.tstep)*Gen.tstep; %[1,TStore]
%Predefine memory for saved pressures
Storage.P = zeros(TStore,Gen.Nx);                               %[TStore,N]
%Store initial pressure
Storage.P(1,:) = State.P;                                       %[TStore,N]

%% - Plotting inputs
Plotting.lwidth_1col = 0.75;
Plotting.Position_1col_matrix = [2.2 1.8 6 4.5];
Plotting.fsize_1col = 7;




end