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
Ny = 1;                                             %[1,1]
N = Nx*Ny;                                          %[1,1]

%Grab pressure
P = reshape(State.P,Nx,Ny);                         %[Nx,1]
%Get density
[Rho,dRhodP] = Density(Flow,P);                     %[Nx,1]
%Get rock transmissibility
[Trans] = Trans_Rock(Flow,Gen);

%Predefine memory
Tx = zeros(Nx+1, Ny);                               %[Nx+1,1]

%% - Transmissibility
%Apply upwind operator to find fluid mobility
%This takes the upwind fluid mobility for the right interface of cell i
Mupx = A.x*Rho/(Flow.muf);                          %[N,1]
%Reshape these 
Mupx = reshape(Mupx, Nx, Ny);                       %[Nx,Ny]

%Update transmissibility to include fluid mobility
%Applies the upwind fluid mobility to the right interface
Tx(2:Nx,:) = Trans.x(2:Nx,:).*Mupx(1:Nx-1,:);       %[Nx+1,Ny]

%Construct matrix 
%Grabs left interface trans (inc mobility) of a cell.
x1 = reshape(Tx(1:Nx,:),N,1);                    	%[N,1]      
%Right interface trans of a cell
x2 = reshape(Tx(2:Nx+1,:),N,1);                     %[N,1] 

% NOTE: Theres a trick here because of the way spdiags cuts off the vectors
% below the diagonal it removes the first cells and above this it removes
% the bottom cells. Because of this we have to swap left and right!
%Collects all transmissibilities
DiagVecs = [-x2,x2+x1,-x1];                   %[N,3]
%Now we end up with a matrix where the row is the cell of interest and the columns are the cells that can either be flowed into or out of
DiagIndx = [-1,0,1];                                %[3,1]
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
dTxdp_L = zeros(Nx+1,Ny);                             %[Nx+1,Ny]
dTxdp_R = zeros(Nx+1,Ny);                             %[Nx+1,Ny]

%Fluid dependence on pressure at right x-interfaces
dMupx = A.x*dRhodP/Flow.muf;                        %[N,1]

%Reshapes the vectors
dMupx = reshape(dMupx, Nx, Ny);                     %[Nx,Ny]

%Multiply the fluid derivative by the appropriate rock transmissibility
dTxdp_L(2:Nx,:) = Trans.x(2:Nx,:).*dMupx(2:Nx,:).*(P(2:Nx,:) - P(1:Nx-1,:));   %[Nx+1,Ny]
dTxdp_R(2:Nx,:) = Trans.x(2:Nx,:).*dMupx(1:Nx-1,:).*(P(2:Nx,:) - P(1:Nx-1,:));   %[Nx+1,Ny]

%Finds the diagonals that will go into the dTransmissibility matrix
%Left
L = reshape(dTxdp_L(1:Nx,:),N,1);                     %[N,1]
%Right
R = -reshape(dTxdp_R(2:Nx+1,:),N,1);                   %[N,1]

%Collect all the dependencies
DiagVecs = [-R, L+R, -L];                       %[N,3]
%Assign to which column they belong
DiagIndx = [-1, 0, 1];                              %[3,1]
%Create a matrix of the dependencies in the correct row
dFTrans = spdiags(DiagVecs,DiagIndx,N,N);           %[N,N]


end

