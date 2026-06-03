%% - Input function
function [Flow, Gen, Plotting, State, Storage] = InputData()
%% - General Parameters
%Final time [sec]
Gen.tf = 10000;                                                 %[1,1]    
%Time step [sec]
Gen.tstep = 100;                                                %[1,1]    
%Tolerance [-]
Gen.tol = 1e-4;                                                 %[1,1]    
%Number of cells in x-direction [-]
Gen.Nx = 5000;                                                  %[1,1]  
%Reservoir length in x-direction [m]
Gen.Lx = 10;                                                    %[1,1]   
%Reservoir length in y-direction [m]
Gen.Ly = 1;                                                     %[1,1]    
%Reservoir height [m]
Gen.Lz = 1;                                                     %[1,1]  
%Fixed boundary pressure on left [Pa]
Gen.PL = 10e6;                                                  %[1,1] 
%Fixed boundary pressure on right [Pa]
Gen.PR = 1e6;                                                   %[1,1] 

%% - Basic Calculations
%Element edge lengths
Gen.dx = Gen.Lx/Gen.Nx;                                         %[1,1]
%Cell locations
Storage.x = linspace(Gen.dx/2,Gen.Lx-Gen.dx/2,Gen.Nx);          %[N,1]

%% - Flow Model
%Permeability [m^2] 
Flow.kx = 1e-12;                                                %[1,1]
%Porosity [-]
Flow.phi = 0.2;                                                 %[1,1]  
%Fluid compressibility [1/Pa]
Flow.cf = 1e-8;                                                 %[1,1]          
%Fluid viscosity [Pa sec]
Flow.muf = .1;                                                  %[1,1]           
%Reference density [kg/m^3]
Flow.Rho0 = 1000;                                               %[1,1]         
%Reference pressure [Pa]
Flow.RhoP = 1e5;                                                %[1,1]  

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