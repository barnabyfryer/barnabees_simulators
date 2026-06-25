function [Q,dQdP] = Add_Wells(Flow,State,Wells)
%% - Well Builder
%Find well source terms and derivatives wrt pressure. All upstreamed.
%Barnaby Fryer - EPFL - 29.03.17

%Calculate perm and derivative
[~,keff,~,~,dkeff_dP] = Perm(Flow,State);                       %[N,1]

%Wellbore flowing pressure [Pa]
Pwf = zeros(length(State.P),1);                                 %[N,1]  
%Only apply locations with wells
Pwf(Wells.Loc_P,1) = Wells.P;                                   %[N,1]

%Wellbore index [m]
WI = zeros(length(State.P),1);                                  %[N,1]  
%Only apply locations with wells
WI(Wells.Loc_P,1) = Wells.WI;                                   %[N,1]

%Upstreamed density of fluid being injected
[Rho(Pwf>=State.P,1),dRhodP(Pwf>=State.P,1)] =...
    Density(Flow,Pwf(Pwf>=State.P,1));                          %[N,1]
dRhodP(Pwf>=State.P,1) = 0;                                     %[N,1]
%Upstreamed density of fluid being produced
[Rho(Pwf<State.P,1),dRhodP(Pwf<State.P,1)] =...
    Density(Flow,State.P(Pwf<State.P,1));                       %[N,1]

%Flow rate [kg/sec]
Q = keff.*(Pwf - State.P).*WI.*Rho/Flow.muf;                    %[N,1]     
%Derivative of flow rate wrt pressure [kg/(Pa*sec)]
dQdP = keff.*WI.*(dRhodP.*(Pwf - State.P) - Rho)/Flow.muf ...
    + dkeff_dP.*(Pwf - State.P).*WI.*Rho/Flow.muf;              %[N,1]

%% - Add fixed rate wells

%Get fixed well data
Q_wells = zeros(length(State.P),1);                             %[N,1]  
Q_wells(Wells.Loc_Q,1) = Wells.Q;                               %[N,1]
%Add this contribution to source terms
Q = Q + Q_wells;                                                %[N,1]

end

