function [A] = Upwind(Gen,State,Trans)
%% - Upwind Calculator
%Returns structure of 2 matrices containg whether a cell (represented by the row)
%is the upstream direction for flow between it and the cell (column) on the right (or top) of row cell
%Barnaby Fryer - EPFL - 29.03.17

P = reshape(State.P,Gen.Nx,Gen.Ny);                         %[Nx,Ny]
Nx = Gen.Nx;                                                %[1,1]
Ny = Gen.Ny;                                                %[1,1]
N = Gen.N;                                                  %[1,1]

%Find velocities at edges. These velocities are still missing all fluid
%properties which need to be upwinded (kr, rho, mu)
U.x = zeros(Nx+1,Ny);                                       %[Nx+1,Ny]
U.y = zeros(Nx,Ny+1);                                       %[Nx,Ny+1]
%Velocity at cell edge in x-direction. If positive flow to right
U.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)).*Trans.x(2:Nx,:);     %[Nx+1,Ny]
%Velocity at cell edge in y-direction. If positive flow to top
U.y(:,2:Ny) = (P(:,1:Ny-1)-P(:,2:Ny)).*Trans.y(:,2:Ny);     %[Nx,Ny+1]

%Use velocity to build upwind operator
%The U.x >0 returns a 1 for all locations with non-negative velocity
R = reshape((U.x(2:Nx+1,:) >= 0), N, 1);                    %[N,1] 
%Then the vector is reshaped to be the size of the reservoir
L = reshape((U.x(1:Nx,:) < 0), N, 1);                       %[N,1]
T = reshape((U.y(:,2:Ny+1) >= 0), N, 1);               	    %[N,1]
B = reshape((U.y(:,1:Ny) < 0), N, 1);                    	%[N,1]

%Stores the vectors containing if the cell is upstream on its right and
%left interfaces
DiagVecs = [R, L];                                          %[N,2]              
%This and the line below build a diagonal matrix. Note!! The vector
%containing if the cell is the upstream direction for the left interface
%is cut-off
DiagIndx = [0,1];                                           %[1,2]                        
%This means that the left directions are shifted such that in the first row
%it will actually be if the second cell is the upstream direction for the
%first cell
A.x = spdiags(DiagVecs,DiagIndx,N,N);                       %[N,1]

%So this matrix will contain whether a cell (represented by the row) is the
%upstream direction for flow between it and the cell (column) on the right 
%(or top) of row cell
DiagVecs = [T, B];                                          %[N,2]   
%This is then done for both x and y directions
DiagIndx = [0,Nx];                                          %[1,2]
A.y = spdiags(DiagVecs,DiagIndx,N,N);                       %[N,N]

end

