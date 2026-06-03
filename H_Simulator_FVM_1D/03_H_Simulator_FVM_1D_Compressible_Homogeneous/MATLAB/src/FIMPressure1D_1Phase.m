function [State] = FIMPressure1D_1Phase(Flow,Gen,State)
%% - Pressure Simulator
%Compressible 2D 1 phase pressure simulator. Uses FIM.
%Barnaby Fryer - EPFL - 29.03.17
    
    %Find the non-fluid rock transmissibility
    [Trans] = Trans_Rock(Flow,Gen);                             %[1,1]
    %Store the state variable from previous time step for use in residual
    State0 = State;
    
    %% - Solve Flow Model
    %Convergence checker
    Converged = 0;                                              %[1,1] 
    %Number of iterations
    it = 1;                                                     %[1,1]
    while Converged == 0
        
        %Build residual
        [Res] =...
            BuildResidual(Flow,Gen,State,State0,Trans);         %[N,1]
        
        %Build Jacobian
        [Jac] = BuildJacobian(Flow,Gen,State,Trans);            %[N,N]
        
        %Solve
        x = -Jac\Res;                                           %[N,1]
        
        %Update Solution
        State.P = State.P + x;                                  %[N,1]
        
        %Build residual
        [Res] =...
            BuildResidual(Flow,Gen,State,State0,Trans);         %[N,1]
        
        %Check Convergence
        %Computes infinite norm of residual
        Norm = norm(Res,inf);                                   %[1,1]                                               
        if Norm < Gen.tol
            %Set to converged
            Converged = 1;                                      %[1,1]           
        else
            %Update iteration tracker
            it = it + 1;                                        %[1,1]                 
        end  
    end

    
    %Find fluid flux
    %Current time step density [kg/m^3]
    [Rho,dRhodP] = Density(Flow,State.P);                           %[N,1]
    %Upwind
    [A] = Upwind(Gen,State,Trans);                                  %[N,N][N,N] 
    %Fluid tansmissivity
    [FTrans,~] = Trans_Fluid(A,dRhodP,Flow,Gen,Rho,Trans);          %[N,N]
    %Find convection into cell
    Conv = FTrans*State.P;                                          %[N,1]
    %Add boundary condition
    [BC,~] = Add_BC(dRhodP,Gen,Rho,State,Trans);                    %[N,1]
    %Fluid flux into cell
    State.flux = -(Conv + BC);                                      %[N,1]

end
