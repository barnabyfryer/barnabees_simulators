%% - Elastic FEM Mechanics Simulator For Quadrilatoral Elements Only 2-D
%Barnaby Fryer 20.04.17

function [State] = M_Simulator_FEM_2D(Gen,Pos,State)


%% - Force vector, determine where to apply forces here
[F] = Build_Force(Gen);

%% - Add pore pressure contribution
[F] = add_pressure(F,Gen,State.P);

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
[stress_e, State.e_vol] = Solve_Stresses(Gen,Pos);                   %[Nn,3]

State.e_vol = reshape(State.e_vol,Gen.Ne,1);
%Effective stress
State.Sig_xx = reshape(stress_e.s_xx,Gen.Ne,1);
State.Sig_yy = reshape(stress_e.s_yy,Gen.Ne,1);
State.Sig_xy = reshape(stress_e.s_xy,Gen.Ne,1);


end
