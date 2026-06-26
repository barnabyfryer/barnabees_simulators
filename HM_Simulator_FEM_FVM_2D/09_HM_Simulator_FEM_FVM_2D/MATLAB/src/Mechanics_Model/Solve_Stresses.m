function [stress_e,e_vol] = Solve_Stresses(Gen,Pos)
%% - FEM code that takes displacements and turns them into stresses
%Barnaby Fryer 27.04.17

%% - Organization
%My nodes are numbered as follows
% 3 - - - 4
% |       |
% |       |
% 1 - - - 2

%And my local directions are
% eta
% ^
% |
% |
% |- - - > xi

%% - Unpack
%Put original node locations in vector form
x = reshape(Pos.x,Gen.Nn,1);                                    %[Nn,1]
y = reshape(Pos.y,Gen.Nn,1);                                    %[Nn,1]
%Horizontal and vertical displacements
u = Pos.u;                                                      %[Nn,1]
v = Pos.v;                                                      %[Nn,1]

%% - Start, decide plane strain or plane stress
AA = (Gen.E/((1+Gen.v)*(1-2*Gen.v)));                           %[1,1]
BB = 1-Gen.v;                                                   %[1,1]
CC = Gen.v;                                                     %[1,1]
DD = ((1-2*Gen.v)/2);                                           %[1,1]

%        _        _
%       | BB CC 0  |
% D = AA| CC BB 0  |
%       |_0  0  DD_|                

%% - Solve for displacements and locations of the nodes of each element
%Predefine memory
E_dx = zeros(Gen.Ne,4);                                         %[Ne,4]
E_dy = zeros(Gen.Ne,4);                                         %[Ne,4]
Loc_x = zeros(Gen.Ne,4);                                        %[Ne,4]
Loc_y = zeros(Gen.Ne,4);                                        %[Ne,4]        

for i = 1:Gen.Ne
    %Row of element
    row = ceil(i/Gen.Nx);                                       %[1,1]
    %Column of element
    col = i - (row-1)*Gen.Nx;                                   %[1,1]
    %Number of element at bottom left
    BL = col + (row-1)*(Gen.Nx+1);                              %[1,1]
    %Number of element at top left
    TL = BL + Gen.Nx + 1;                                       %[1,1]
    
    %% - Displacements of the nodes of each element
    %x-displacement
    %x-displacement of bottom left node
    E_dx(i,1) = u(BL,1);                                        %[Ne,4]
    %x-displacement of bottom right node
    E_dx(i,2) = u(BL+1,1);                                      %[Ne,4]
    %x-displacement of top left node
    E_dx(i,3) = u(TL,1);                                        %[Ne,4]
    %x-displacement of top right node
    E_dx(i,4) = u(TL+1,1);                                      %[Ne,4]
    
    %y-displacement
    %y-displacement of bottom left node
    E_dy(i,1) = v(BL,1);                                        %[Ne,4]
    %y-displacement of bottom right node
    E_dy(i,2) = v(BL+1,1);                                      %[Ne,4]
    %y-displacement of top left node
    E_dy(i,3) = v(TL,1);                                        %[Ne,4]
    %y-displacement of top right node
    E_dy(i,4) = v(TL+1,1);                                      %[Ne,4]
    
    %% - Original global coordinates of each node of each element
    %x-location
    %x-location of bottom left node
    Loc_x(i,1) = x(BL,1);                                       %[Ne,4]
    %x-location of bottom right node
    Loc_x(i,2) = x(BL+1,1);                                     %[Ne,4]
    %x-location of top left node
    Loc_x(i,3) = x(TL,1);                                       %[Ne,4]
    %x-location of top right node
    Loc_x(i,4) = x(TL+1,1);                                     %[Ne,4]
    
    %y-location
    %y-location of bottom left node
    Loc_y(i,1) = y(BL,1);                                       %[Ne,4]
    %y-location of bottom right node
    Loc_y(i,2) = y(BL+1,1);                                     %[Ne,4]
    %y-location of top left node
    Loc_y(i,3) = y(TL,1);                                       %[Ne,4]
    %y-location of top right node
    Loc_y(i,4) = y(TL+1,1);                                     %[Ne,4]
end

%% - D Matrix
D = zeros(3,3);                                                 %[4,4]
D(1,1) = BB;                                                    %[4,4]
D(1,2) = CC;                                                    %[4,4]
D(2,1) = CC;                                                    %[4,4]
D(2,2) = BB;                                                    %[4,4]
D(3,3) = DD;                                                    %[4,4]

D = D*AA;                                                       %[4,4]

%% - Solve for strains and stresses 

%Local coordinate which corresponds to y direction
eta = 0;                                                        %[1,1]

Strain = zeros(Gen.Ne,3);                                       %[Ne,3]
Sigma = zeros(Gen.Ne,3);                                        %[Ne,3]

%Loop over all elements
for i = 1:Gen.Ne
    
    %Get the local coordinates of the nodes corresponding to this element
    coord(:,1) = Loc_x(i,:)';                                   %[4,2]
    coord(:,2) = Loc_y(i,:)';                                   %[4,2]
    
    %Find the displacements of nodes corresponding to this element
    disp = zeros(1,8);                                          %[1,8]
    disp(:,1:2:7) = E_dx(i,:);                                  %[1,8]
    disp(:,2:2:8) = E_dy(i,:);                                  %[1,8]
    
    %Loop over local coordinate xi (which corresponds to x direction)
    % for xi = -1:.05:1
    for xi = 0
        N = zeros(4,1);                                         %[4,1]
        dN = zeros(2,4);                                        %[2,4]
        B = zeros(3,8);                                         %[3,8]
            
        %% - Value of shape functions
        N(1,1) = 1/4*(1-xi)*(1-eta);                            %[4,1]                    
        N(2,1) = 1/4*(1+xi)*(1-eta);                            %[4,1]  
        N(3,1) = 1/4*(1-xi)*(1+eta);                            %[4,1]  
        N(4,1) = 1/4*(1+xi)*(1+eta);                            %[4,1]  
        
        %% - Value of derivatives of shape functions
        %Written as dN(Direction,Node)
        %Node 1
        dN(1,1) = -1/4*(1-eta);                                 %[2,4]
        dN(2,1) = -1/4*(1-xi);                                  %[2,4]
        %Node 2
        dN(1,2) = 1/4*(1-eta);                                  %[2,4]
        dN(2,2) = -1/4*(1+xi);                                  %[2,4]
        %Node 3
        dN(1,3) = -1/4*(1+eta);                                 %[2,4]
        dN(2,3) = 1/4*(1-xi);                                   %[2,4]
        %Node 4
        dN(1,4) = 1/4*(1+eta);                                  %[2,4]
        dN(2,4) = 1/4*(1+xi);                                   %[2,4]
        
        %% - Jacobian calculations
        %Solve for the Jacobian
        Jac = dN*coord;                                         %[2,2]
        %Find the determinant of the Jacobian
        DJ = det(Jac);                                          %[1,1]
        %Find derivatives for B matrix
        deriv = Jac\dN;                                         %[2,4]
        
        %% - Formulate B Matrix
        %dN1/dx
        B(1,1) = deriv(1,1);                                    %[3,8]
        %dN1/dy
        B(3,1) = deriv(2,1);                                    %[3,8]
        %dN1/dy
        B(2,2) = deriv(2,1);                                    %[3,8]
        %dN1/dx
        B(3,2) = deriv(1,1);                                    %[3,8]
        %dN2/dx
        B(1,3) = deriv(1,2);                                    %[3,8]
        %dN2/dy
        B(3,3) = deriv(2,2);                                    %[3,8]
        %dN2/dy
        B(2,4) = deriv(2,2);                                    %[3,8]
        %dN2/dx
        B(3,4) = deriv(1,2);                                    %[3,8]
        %dN3/dx
        B(1,5) = deriv(1,3);                                    %[3,8]
        %dN3/dy
        B(3,5) = deriv(2,3);                                    %[3,8]
        %dN3/dy
        B(2,6) = deriv(2,3);                                    %[3,8]
        %dN3/dx
        B(3,6) = deriv(1,3);                                    %[3,8]
        %dN4/dx
        B(1,7) = deriv(1,4);                                    %[3,8]
        %dN4/dy
        B(3,7) = deriv(2,4);                                    %[3,8]
        %dN4/dy
        B(2,8) = deriv(2,4);                                    %[3,8]
        %dN4/dx
        B(3,8) = deriv(1,4);                                    %[3,8]
        
        %% - Strain and Stress calculations
        %Calculate the strains of the element, positive in compression
        Strain(i,:) = -(B*transpose(disp(1,:)))';               %[Ne,3]
        %Calculate the effective stresses of the element
        Sigma(i,:) = (D*Strain(i,:)')';                         %[Ne,3]
    end
    
    %Find volumetric strain
    e_vol = reshape(Strain(:,1),Gen.Nx,Gen.Ny) + reshape(Strain(:,2),Gen.Nx,Gen.Ny);
	
	%Effective stress changes
    stress_e.s_xx = reshape(Sigma(:,1),Gen.Nx,Gen.Ny);
    stress_e.s_yy = reshape(Sigma(:,2),Gen.Nx,Gen.Ny);
    stress_e.s_xy = reshape(Sigma(:,3),Gen.Nx,Gen.Ny);

end

end

