function [k_global] = Global_Stiffness_Builder(Gen)

%Predefine memory of global stiffness matrix
k_global = zeros(2*Gen.Nn,2*Gen.Nn);                        %[2*Nn,2*Nn]

%Build global stiffness matrix
for i = 1:Gen.Ne
    %Find nodes in element
    Nodes = find(Gen.Ref(i,:),4);                           %[1,4]
    %Prepare local stiffness matrix of element
    %Plane strain
    [kl] = Local_Stiffness_PStrain(Gen);              %[8,8]
    %Plane stress
    %[kl] = Local_Stiffness_PStress(Gen);              %[8,8] 
    
    %Add element stiffness matrix to global stiffness matrix
    %Apply local stiffness matrix contribution to each node of the element
    for j = 1:8
        row = Nodes(floor((j-1)/2) + 1);                    %[1,1]
        loc = 2*(row) - 1 + rem(j-1,2);                     %[1,1]
        
        %Pick all even and odd columns of row of k_global
        odd = kl(j,1:2:length(kl(j,:)));                    %[1,4]
        even = kl(j,2:2:length(kl(j,:)));                   %[1,4]
        
        %Find the columns they correspond to
        oddI = Nodes*2 - 1;                                 %[1,4]
        evenI = Nodes*2;                                    %[1,4]
        
        %Add local stiffness matrix to global stiffness matrix
        k_global(loc,oddI) = k_global(loc,oddI) + odd;      %[2*Nn,2*Nn]
        k_global(loc,evenI) = k_global(loc,evenI) + even;   %[2*Nn,2*Nn]
        
    end 
end

end