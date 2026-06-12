function [Jac] = BuildJacobian(Flow,Gen,State,Trans,Wells)
%% - Jacobian Builder
%Builds FIM Jacobian
%Barnaby Fryer - EPFL - 29.03.17

%Basic grid parameters
Nx = Gen.Nx;                                            %[1,1]
Ny = 1;                                                 %[1,1]
N = Nx*Ny;                                              %[1,1]
dx = Gen.Lx/Nx;                                         %[1,1]
dy = Gen.Ly/Ny;                                         %[1,1]
V = dx*dy*Gen.Lz;                                       %[1,1]

%Returns density and derivative of density wrt pressure of each cell
[Rho,dRhodP] = Density(Flow,State.P);                   %[N,1] [N,1]

%Calculate porosity and derivative
[phi,dphidp] = PhiCalc(Flow,State);                	    %[N,1]

 %Returns structure of 2 matrices containg whether a cell (represented by the row) is the upstream direction for flow between it and the cell (column) on the right (or top) of row cell
[A] = Upwind(Gen,State,Trans);                          %2x[N,N] 

%Returns the transmissibilities between the two cells
[FTrans,dFTrans] = Trans_Fluid(A,Flow,Gen,State);
%The sum of all transmissibilities is also located in the main diagonal
%basically this is so you can multiply FTrans*P (where P is a vector of
%pressure) and get the equations T_R * (P-P_R) + T_L * (P-P_L) etc
%It also returns the derivative of each cell's transmissibility wrt to the
%pressure of every cell, multiplied by the difference in pressure of those
%two cells. This is effectively the effect compressibility has

%% - Derivative of accumulation terms
Acc = (dRhodP.*phi + Rho.*dphidp)*V/Gen.tstep;          %[N,1]

%% - Derivative of convection terms
%Due to change in pressure
Conv = FTrans;                                          %[N,N]

%Due to compressibility and permeability changes
Comp = dFTrans;                                         %[N,N]

%% - Derivative of source terms
[~,dQdP] = Add_Wells(Flow,State,Wells);                 %[N,1]

%% - Combine terms
%Adds main diagonal terms from accumulation and wells
DiagVecs = Acc - dQdP;                                  %[N,1]         
%Chooses main diagonal
DiagIndx = 0;                                           %[1,1]
%Adds derivatives of convection terms
Jac = spdiags(DiagVecs,DiagIndx,N,N) + Conv + Comp;     %[N,N] 

end

