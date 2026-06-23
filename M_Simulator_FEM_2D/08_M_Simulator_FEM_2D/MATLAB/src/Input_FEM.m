function [Gen, Plotting, Pos] = Input_FEM()


%% - General Inputs
%Elements in x direction [-]
Gen.Nx = 21;                                         	%[1,1]
%Elements in y direction [-]
Gen.Ny = 21;                                        	%[1,1]
%Reservoir lengths [m]
Gen.Lx = 30;                                        	%[1,1]
Gen.Ly = 10;                                            %[1,1]
%Young's Modulus [Pa]
Gen.E = 1e9;                                            %[1,1]
%Poisson's Ratio [-]
Gen.v = .3;                                             %[1,1]

%% - Basic Calculations
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

%% - Plotting inputs
Plotting.lwidth_1col = 0.75;
Plotting.Position_1col_matrix = [2.2 1.8 6 4.5];
Plotting.fsize_1col = 7;

end

