function [FTrans,dFTrans] = Trans_Fluid(A,Flow,Gen,State)
%% - Full Transmissibility Builder
%Returns the transmissibilities between the two cells
%The sum of all transmissibilities is also located in the main diagonal
%basically this is so you can multiply FTrans*P (where P is a vector of
%pressure) and get the equations T_R * (P-P_R) + T_L * (P-P_L) etc
%It also returns the derivative of each cell's transmissibility wrt to the
%pressure of every cell, multiplied by the difference in pressure of those
%two cells. This is effectively the effect compressibility has.
%Barnaby Fryer - EPFL - 29.03.17

%Basic grid properties
Nx = Gen.Nx;                                        %[1,1]
Ny = Gen.Ny;                                        %[1,1]
N = Gen.N;                                          %[1,1]

%Grab pressure
P = reshape(State.P,Nx,Ny);                         %[Nx,1]
%Get density
[Rho,dRhodP] = Density(Flow,State.P);               %[Nx,1]
%Get rock transmissibility
[Trans,dTransdP_LB,dTransdP_RT] = Trans_Rock(Flow,Gen,State);

%Predefine memory
Tx = zeros(Nx+1, Ny);                               %[Nx+1,Ny]
Ty = zeros(Nx, Ny+1);                               %[Nx,Ny+1]

%% - Transmissibility
%Apply upwind operator to find fluid mobility
%This takes the upwind fluid mobility for the right interface of cell i
Mupx = A.x*Rho/(Flow.muf);                          %[N,1]
Mupy = A.y*Rho/(Flow.muf);                          %[N,1]
%Reshape these 
Mupx = reshape(Mupx, Nx, Ny);                       %[Nx,Ny]
Mupy = reshape(Mupy, Nx, Ny);                       %[Nx,Ny]

%Update transmissibility to include fluid mobility
%Applies the upwind fluid mobility to the right interface
Tx(2:Nx,:) = Trans.x(2:Nx,:).*Mupx(1:Nx-1,:);       %[Nx+1,Ny]
%Does the same thing for top interface
Ty(:,2:Ny) = Trans.y(:,2:Ny).*Mupy(:,1:Ny-1);       %[Nx,Ny+1]

%Construct matrix 
%Grabs left interface trans (inc mobility) of a cell.
x1 = reshape(Tx(1:Nx,:),N,1);                    	%[N,1]      
%Right interface trans of a cell
x2 = reshape(Tx(2:Nx+1,:),N,1);                     %[N,1] 
%Bottom
y1 = reshape(Ty(:,1:Ny),N,1);                       %[N,1]
%Top
y2 = reshape(Ty(:,2:Ny+1),N,1);                     %[N,1]

% NOTE: Theres a trick here because of the way spdiags cuts off the vectors
% below the diagonal it removes the first cells and above this it removes
% the bottom cells. Because of this we have to swap left and right!
%Collects all transmissibilities
DiagVecs = [-y2,-x2,y2+x2+x1+y1,-x1,-y1];           %[N,5]
%Now we end up with a matrix where the row is the cell of interest and the columns are the cells that can either be flowed into or out of
DiagIndx = [-Nx,-1,0,1,Nx];                         %[5,1]
%The values are the transmissibilities between the two cells
FTrans = spdiags(DiagVecs,DiagIndx,N,N);            %[N,N]

%The sum of all transmissibilities is also located in the main diagonal
%basically this is so you can multiply FTrans*P (where P is a vector of
%pressure) and get the equations T_R * (P-P_R) + T_L * (P-P_L) etc

%% - dTransdP, part 1 (density dependence on pressure)
%The derivative of each cell's transmissibility wrt to the
%pressure of every cell, multiplied by the difference in pressure of those
%two cells. This is effectively the effect compressibility has.

%Predefines memory for derivative of trans wrt pressure
%Derivative due to dependence of x interface on pressure
dTxdp_L = zeros(Nx+1,Ny);                               %[Nx+1,Ny]
dTxdp_R = zeros(Nx+1,Ny);                               %[Nx+1,Ny]
%And for y-direction interfaces
dTydp_B = zeros(Nx,Ny+1);                               %[Nx,Ny+1]
dTydp_T = zeros(Nx,Ny+1);                         	    %[Nx,Ny+1]

%Fluid dependence on pressure at right x-interfaces
dMupx = A.x*dRhodP/Flow.muf;                            %[N,1]
dMupy = A.y*dRhodP/Flow.muf;                            %[N,1]

%Reshapes the vectors
dMupx = reshape(dMupx, Nx, Ny);                         %[Nx,Ny]
dMupy = reshape(dMupy, Nx, Ny);                         %[Nx,Ny]

%Multiply the fluid derivative by the appropriate rock transmissibility
dTxdp_L(2:Nx,:) = Trans.x(2:Nx,:).*dMupx(2:Nx,:).*(P(2:Nx,:) - P(1:Nx-1,:));        %[Nx+1,Ny]
dTxdp_R(2:Nx,:) = Trans.x(2:Nx,:).*dMupx(1:Nx-1,:).*(P(2:Nx,:) - P(1:Nx-1,:));      %[Nx+1,Ny]
dTydp_B(:,2:Ny) = Trans.y(:,2:Ny).*dMupy(:,2:Ny).*(P(:,2:Ny) - P(:,1:Ny-1));        %[Nx,Ny+1]
dTydp_T(:,2:Ny) = Trans.y(:,2:Ny).*dMupy(:,1:Ny-1).*(P(:,2:Ny) - P(:,1:Ny-1));      %[Nx,Ny+1]

%Finds the diagonals that will go into the dTransmissibility matrix
%Left
L = reshape(dTxdp_L(1:Nx,:),N,1);                   %[N,1]
%Right
R = -reshape(dTxdp_R(2:Nx+1,:),N,1);                %[N,1]
%Bottom
B = reshape(dTydp_B(:,1:Ny),N,1);             	    %[N,1]
%Top
T = -reshape(dTydp_T(:,2:Ny+1),N,1);         	    %[N,1]

%Collect all the dependencies
DiagVecs = [-T, -R, T+L+R+B, -L, -B];               %[N,5]
%Assign to which column they belong
DiagIndx = [-Nx, -1, 0, 1, Nx];                     %[5,1]
%Create a matrix of the dependencies in the correct row
dFTrans_1 = spdiags(DiagVecs,DiagIndx,N,N);         %[N,N]

%% - dTransdP, part 2 (permeability dependence on pressure)
%Predefines memory for derivative of trans wrt pressure
%Derivative due to dependence of x interface on pressure
dTxdp_L = zeros(Nx+1,Ny);                           %[Nx+1,Ny]
dTxdp_R = zeros(Nx+1,Ny);                           %[Nx+1,Ny]
%For derivatives in the y direction
dTydp_B = zeros(Nx,Ny+1);                           %[Nx,Ny+1]
dTydp_T = zeros(Nx,Ny+1);                           %[Nx,Ny+1]

%Fluid dependence on pressure at right x-interfaces
Mupx = A.x*Rho/Flow.muf;                            %[N,1]
Mupy = A.y*Rho/Flow.muf;                            %[N,1]

%Reshapes the vectors
Mupx = reshape(Mupx, Nx, Ny);                       %[Nx,Ny]
Mupy = reshape(Mupy, Nx, Ny);                       %[Nx,Ny]

%Multiply the fluid derivative by the appropriate rock transmissibility
dTxdp_L(2:Nx,:) = dTransdP_LB.x(2:Nx,:).*Mupx(2:Nx,:).*(P(2:Nx,:) - P(1:Nx-1,:));       %[Nx+1,Ny]
dTxdp_R(2:Nx,:) = dTransdP_RT.x(2:Nx,:).*Mupx(1:Nx-1,:).*(P(2:Nx,:) - P(1:Nx-1,:));     %[Nx+1,Ny]
dTydp_B(:,2:Ny) = dTransdP_LB.y(:,2:Ny).*Mupy(:,2:Ny).*(P(:,2:Ny) - P(:,1:Ny-1));       %[Nx,Ny+1]
dTydp_T(:,2:Ny) = dTransdP_LB.y(:,2:Ny).*Mupy(:,1:Ny-1).*(P(:,2:Ny) - P(:,1:Ny-1));     %[Nx,Ny+1]

%Finds the diagonals that will go into the dTransmissibility matrix
%Left
L = reshape(dTxdp_L(1:Nx,:),N,1);                   %[N,1]
%Right
R = -reshape(dTxdp_R(2:Nx+1,:),N,1);                %[N,1]
%Bottom
B = reshape(dTydp_B(:,1:Ny),N,1);                   %[N,1]
%Top
T = -reshape(dTydp_T(:,2:Ny+1),N,1);                %[N,1]

%Collect all the dependencies
DiagVecs = [-T, -R, T+L+R+B, -L -B];                %[N,5]
%Assign to which column they belong
DiagIndx = [-Nx, -1, 0, 1, Nx];                     %[5,1]
%Create a matrix of the dependencies in the correct row
dFTrans = spdiags(DiagVecs,DiagIndx,N,N) + dFTrans_1;           %[N,N]

end

