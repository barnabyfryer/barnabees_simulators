function [F] = Build_Force(Gen)

%Where to apply the force (on all top nodes)
A = (2*Gen.Nn-Gen.Nx*2+2):2:2*Gen.Nn-2;       
%On all right nodes
B = 4*(Gen.Nx+1)-1:2*(Gen.Nx+1):2*Gen.Nn-2*(Gen.Nx+1);
%On all bottom nodes
C = 4:2:2*Gen.Nx;
%On all left nodes
D = 2*(Gen.Nx+1)+1:2*(Gen.Nx+1):2*Gen.Nn-4*(Gen.Nx+1)+1;
%Set up force vector
F = zeros(2*Gen.Nn,1);                                      %[2*Nn,1]

%Apply force (per unit of out of plane thickness) [N/m]
fy = -0e7;
fx = -0e7;

%On top nodes
F(A,1) = fy;
%On right nodes
F(B,1) = fx;
%On bottom nodes
% F(C,1) = -fy;
%On left nodes
% F(D,1) = -fx;

%End of force 
%On top left node y direction
F(2*Gen.Nn-Gen.Nx*2,1) = fy/2;
%On top left node x direction
% F(2*Gen.Nn-Gen.Nx*2-1,1) = -fx/2;
%On top right node y direction
F(2*Gen.Nn,1) = fy/2;
%On top right node x direction
F(2*Gen.Nn-1,1) = fx/2;
%On bottom right node y dirction
% F(2*(Gen.Nx+1),1) = -fy/2;
%On bottom right node x direction
F(2*(Gen.Nx+1)-1,1) = fx/2;
%On bottom left node y direction
% F(2,1) = -fy/2;
%On bottom left node x direction
% F(1,1) = -fx/2;

%Add pressure contributions to force vector
% [F] = Add_Pressure(F,Gen,Res.P);                            %[2*Nn,1]

end