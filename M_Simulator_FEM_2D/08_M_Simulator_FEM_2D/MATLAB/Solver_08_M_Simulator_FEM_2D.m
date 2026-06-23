%% - Elastic FEM Mechanics Simulator For Quadrilatoral Elements Only 2-D
%Barnaby Fryer 20.04.17

clc
clear all
close all

%% - Input Data
[Gen, Plotting, Pos] = Input_FEM();

%% - Force vector, determine where to apply forces here
[F] = Build_Force(Gen);

%% - Stiffness Matrix
[k_global] = Global_Stiffness_Builder(Gen);

%% - Boundary Conditions
%Apply fixes to global matrix
[F,k_global] = Fix_x(F,k_global,Gen.nodesx);
[F,k_global] = Fix_y(F,k_global,Gen.nodesy);

%% - Solve
Pos.dx = k_global\F;                                            %[2*Nn,1]

%% - Process displacements
%Displacement in x-direction
Pos.u = Pos.dx(1:2:length(Pos.dx),1);                           %[Nn,1]
%Displacement in y-direction
Pos.v = Pos.dx(2:2:length(Pos.dx),1);                           %[Nn,1]

%Displacement in x-direction reshaped
Pos.du = reshape(Pos.u,Gen.Nx+1,Gen.Ny+1);                      %[Nx+1,Ny+1]
%Displacement in y-direction reshaped
Pos.dv = reshape(Pos.v,Gen.Nx+1,Gen.Ny+1);                      %[Nx+1,Ny+1]

%New node locations
Pos.x_new = Pos.x + Pos.du;                                     %[Nx+1,Ny+1]
Pos.y_new = Pos.y + Pos.dv;                                     %[Nx+1,Ny+1]

%% - Solve for stresses
%Solve for stresses
[Sig] = Solve_Stresses(Gen,Pos);                                %[Nn,3]

%% - Plot
Plotting_sim(F,Gen,Plotting,Pos,Sig);
