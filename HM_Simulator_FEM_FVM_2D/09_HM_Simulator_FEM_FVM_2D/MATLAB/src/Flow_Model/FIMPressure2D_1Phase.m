function [State] = FIMPressure2D_1Phase(Flow,Gen,State,Wells)
%% - Pressure Simulator
%Compressible 2D 1 phase pressure simulator. Uses FIM.
%Barnaby Fryer - EPFL - 29.03.17
    
    
    %Store the state variable from previous time step for use in residual
    State0 = State;
    
    %% - Solve Flow Model
    %Convergence checker
    Converged = 0;                                              %[1,1] 
    %Number of iterations
    it = 1;                                                     %[1,1]
    while Converged == 0

        %Find the non-fluid rock transmissibility
        [Trans,~,~] = Trans_Rock(Flow,Gen,State);               %2x[Nx,1]
        
        %Build residual
        [Res] =...
            BuildResidual(Flow,Gen,State,State0,Trans,Wells);   %[N,1]
        
        %Build Jacobian
        [Jac] = BuildJacobian(Flow,Gen,State,Trans,Wells);      %[N,N]
        
        %Solve
        x = -Jac\Res;                                           %[N,1]
        
        %Update Solution
        State.P = State.P + x;                                  %[N,1]
        
        %Build residual
        [Res] =...
            BuildResidual(Flow,Gen,State,State0,Trans,Wells);   %[N,1]
        
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
    %Upwind
    [A] = Upwind(Gen,State,Trans);                              %[N,N][N,N] 
    %Fluid tansmissivity
    [FTrans,~] = Trans_Fluid(A,Flow,Gen,State);                 %[N,N]
    %Find source terms
    [Q,~] = Add_Wells(Flow,State,Wells);                        %[N,1]
    %Find convection into cell
    State.flux = -FTrans*State.P + Q;                           %[N,1]
    %Calculate perm and derivative
    [State,~] = Perm(Flow,State);                	            %[N,1]
    %Calculate porosity and derivative
    [State.phi,~] = PhiCalc(Flow,State);                	    %[N,1]

end
