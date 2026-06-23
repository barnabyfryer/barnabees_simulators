function [k_local] = Local_Stiffness_PStrain(Gen)

%Length of local coordinate
a = Gen.dx;                                 %[1,1]
b = Gen.dy;                                 %[1,1]

%Poisson's ratio
v = Gen.v;                                  %[1,1]

%Define constants specific to plane stress
A = (Gen.E/((1+v)*(1-2*v)));                %[1,1]
B = 1-v;                                    %[1,1]
C = v;                                      %[1,1]
D = ((1-2*v)/2);                            %[1,1]

%Predefine memory
k_local = zeros(8,8);                       %[8,8]

%% - Build each stiffness matrix component
%Row 1
k_local(1,1) = B*b/(3*a) + D*(a/(3*b));
k_local(1,2) = C*(1/4) + D*(1/4);
k_local(1,3) = -B*b/(3*a) + D*(a/(6*b));
k_local(1,4) = C*(1/4) + D*(-1/4);
k_local(1,5) = B*b/(6*a) + D*(-a/(3*b));
k_local(1,6) = C*(-1/4) + D*(1/4);
k_local(1,7) = -B*b/(6*a) + D*(-a/(6*b));
k_local(1,8) = C*(-1/4) + D*(-1/4);

%Row 2
k_local(2,1) = C*(1/4) + D*(1/4);
k_local(2,2) = B*a/(3*b) + D*(b/(3*a));
k_local(2,3) = C*(-1/4) + D*(1/4);
k_local(2,4) = B*a/(6*b) + D*(-b/(3*a));
k_local(2,5) = C*(1/4) + D*(-1/4);
k_local(2,6) = -B*a/(3*b) + D*(b/(6*a));
k_local(2,7) = C*(-1/4) + D*(-1/4);
k_local(2,8) = -B*a/(6*b) + D*(-b/(6*a));

%Row 3
k_local(3,1) = -B*b/(3*a) + D*(a/(6*b));
k_local(3,2) = C*(-1/4) + D*(1/4);
k_local(3,3) = B*b/(3*a) + D*(a/(3*b));
k_local(3,4) = C*(-1/4) + D*(-1/4);
k_local(3,5) = -B*b/(6*a) + D*(-a/(6*b));
k_local(3,6) = C*(1/4) + D*(1/4);
k_local(3,7) = B*b/(6*a) + D*(-a/(3*b));
k_local(3,8) = C*(1/4) + D*(-1/4);

%Row 4
k_local(4,1) = C*(1/4) + D*(-1/4);
k_local(4,2) = B*a/(6*b) + D*(-b/(3*a));
k_local(4,3) = C*(-1/4) + D*(-1/4);
k_local(4,4) = B*a/(3*b) + D*(b/(3*a));
k_local(4,5) = C*(1/4) + D*(1/4);
k_local(4,6) = -B*a/(6*b) + D*(-b/(6*a));
k_local(4,7) = C*(-1/4) + D*(1/4);
k_local(4,8) = -B*a/(3*b) + D*(b/(6*a));

%Row 5
k_local(5,1) = B*b/(6*a) + D*(-a/(3*b));
k_local(5,2) = C*(1/4) + D*(-1/4);
k_local(5,3) = -B*b/(6*a) + D*(-a/(6*b));
k_local(5,4) = C*(1/4) + D*(1/4);
k_local(5,5) = B*b/(3*a) + D*(a/(3*b));
k_local(5,6) = C*(-1/4) + D*(-1/4);
k_local(5,7) = -B*b/(3*a) + D*(a/(6*b));
k_local(5,8) = C*(-1/4) + D*(1/4);

%Row 6
k_local(6,1) = C*(-1/4) + D*(1/4);
k_local(6,2) = -B*a/(3*b) + D*(b/(6*a));
k_local(6,3) = C*(1/4) + D*(1/4);
k_local(6,4) = -B*a/(6*b) + D*(-b/(6*a));
k_local(6,5) = C*(-1/4) + D*(-1/4);
k_local(6,6) = B*a/(3*b) + D*(b/(3*a));
k_local(6,7) = C*(1/4) + D*(-1/4);
k_local(6,8) = B*a/(6*b) + D*(-b/(3*a));

%Row 7
k_local(7,1) = -B*b/(6*a) + D*(-a/(6*b));
k_local(7,2) = C*(-1/4) + D*(-1/4);
k_local(7,3) = B*b/(6*a) + D*(-a/(3*b));
k_local(7,4) = C*(-1/4) + D*(1/4);
k_local(7,5) = -B*b/(3*a) + D*(a/(6*b));
k_local(7,6) = C*(1/4) + D*(-1/4);
k_local(7,7) = B*b/(3*a) + D*(a/(3*b));
k_local(7,8) = C*(1/4) + D*(1/4);

%Row 8
k_local(8,1) = C*(-1/4) + D*(-1/4);
k_local(8,2) = -B*a/(6*b) + D*(-b/(6*a));
k_local(8,3) = C*(1/4) + D*(-1/4);
k_local(8,4) = -B*a/(3*b) + D*(b/(6*a));
k_local(8,5) = C*(-1/4) + D*(1/4);
k_local(8,6) = B*a/(6*b) + D*(-b/(3*a));
k_local(8,7) = C*(1/4) + D*(1/4);
k_local(8,8) = B*a/(3*b) + D*(b/(3*a));

%Apply A
k_local = k_local*A;

end


