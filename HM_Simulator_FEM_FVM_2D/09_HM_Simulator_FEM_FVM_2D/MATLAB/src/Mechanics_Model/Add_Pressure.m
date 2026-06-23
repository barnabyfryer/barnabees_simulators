function [F] = Add_Pressure(F,Gen,P)

P = reshape(P,Gen.Nn,1);

for i = 1:Gen.Ne
    %Find nodes in element
    Nodes = find(Gen.Ref(i,:),4);                           %[1,4]
    
    
    P_avg = sum(P(Nodes))/4;
    
    F_app = P_avg*(Gen.dx*Gen.dy);
    
    
    NodesB = Nodes(1,1:2);
    NodesT = Nodes(1,3:4);
    NodesL = Nodes(1,1:2:3);
    NodesR = Nodes(1,2:2:4);
    
    %X-direction forces
    %Positive forces
    F(NodesR*2-1,1) = F(NodesR*2-1,1) + F_app;
    %Negative forces
    F(NodesL*2-1,1) = F(NodesL*2-1,1) - F_app;
    
    %Y-direction forces
    %Positive forces
    F(NodesT*2,1) = F(NodesT*2,1) + F_app;
    %Negative forces
    F(NodesB*2,1) = F(NodesB*2,1) - F_app;
    
    
    
end

end

