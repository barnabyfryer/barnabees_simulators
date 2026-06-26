function [Res] = BuildResidual(Flow,Gen,State,State0,Trans,Wells)
%% - Residual Builder
%Builds FIM residual
%Barnaby Fryer - EPFL - 29.03.17

%% - Basic Calculations
%Returns structure of 2 matrices containg whether a cell (represented by
%the row) is the upstream direction for flow between it and the cell 
%(column) on the right (or top) of row cell      
[A] = Upwind(Gen,State,Trans);                                  %[N,N][N,N] 

%Old time step density [kg/m^3]
[RhoOld,~] = Density(Flow,State0.P);                            %[N,1]
%Current iteration density [kg/m^3]
[Rho,~] = Density(Flow,State.P);                                %[N,1]

%Calculate porosity of last time step
[phiOld,~] = PhiCalc(Flow,State0);                	            %[N,1]
%Calculate porosity and derivative
[phi,~] = PhiCalc(Flow,State);                	                %[N,1]

%Transmissibility matrix
[FTrans,~] = Trans_Fluid(A,Flow,Gen,State);                     %[N,N]

%Length of cell in y-direction [m]
dy = Gen.Ly/1;                                                  %[1,1]                
%Length of cell in y-direction [m]
dx = Gen.Lx/Gen.Nx;                                             %[1,1]  
%Cross-sectional area in x-direction [m^2]
Ax = dy*Gen.Lz;                                                 %[1,1]    
%Cell volume [m^3]
V = Ax*dx;                                                      %[1,1]

%% - Accumulation Term [kg/sec]
Acc = (Rho.*phi - RhoOld.*phiOld)*(V/Gen.tstep);                %[N,1]

%% - Convection Term [kg/sec]
Conv = FTrans*State.P;                                          %[N,1]

%% - Source Term [kg/sec]
[Q,~] = Add_Wells(Flow,State,Wells);                            %[N,1]

%% - Find Residual [kg/sec]
Res = Acc + Conv - Q;                                           %[N,1]

end

