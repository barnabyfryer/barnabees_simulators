clc
clear all
close all

%% - About
%This reservoir simulator uses a FVM formulation to solve the continuity of
%mass balance equation in 1-D. The boundary conditions are no flow
%at the edges (Neumann). It uses a pressure-dependent, heterogeneous permeability and porosity.
%There is no gravity and the simulator is single phase. The fluid is
%considered to be slightly compressible. 
%That flow model is coupled to a FEM mechanics model which is elastic and assumes small strains.
%The coupling is two-way and iterative.

%Barnaby Fryer, 2026

%% - Inputs
[Flow, Gen, Plotting, State, Storage, Wells] = InputFlowData();
[Gen, Pos, State, Storage] = Input_FEM(Flow, Gen, State, Storage);

%% - Run Simulation
while State.t < Gen.tf

    %Store current iteration
    it = 0;
    %Initialize error
    err = 1;
    %Store the previous time step for use in residual
    State0 = State;
    %Store the previous iteration for density update after mechanical model
    State_phi = State;
    while err > Gen.tol_all
        %Save state prior to iteration
        P0 = State.P;
        e_vol0 = State.e_vol;

        %% - Solve For Pressure
        [State] = FIMPressure2D_1Phase(Flow,Gen,State,State0,State_phi,Wells);

        %Store the previous iteration for density update after mechanical model
        State_phi = State;

        %% - Solve for stresses and volumetric strain
        [State] = M_Simulator_FEM_2D(Gen,Pos,State);

        %% - Check error
        dP = State.P - P0;
        State.errP = abs(dP) / max(max(abs(State.P)),1e5);
        State.erre = abs((State.e_vol - e_vol0))./max(max(abs(State.e_vol)),1e-12);
        %Find error
        err = max([max(State.errP), max(State.erre)]);
        %Update iteration
        it = it + 1;

    end

    %Store results
    if ~isempty(find(State.t == Storage.TStorage,1))
        State.step = State.step + 1;
        Storage.P(State.step,:) = State.P;
        Storage.flux(State.step,:) = State.flux;
        Storage.kx(State.step,:) = State.kx;
        Storage.ky(State.step,:) = State.ky;
        Storage.phi(State.step,:) = State.phi;
        Storage.e_vol(State.step,:) = State.e_vol;
        Storage.errP(State.step,:) = State.errP;
        Storage.erre(State.step,:) = State.erre;
    end

    State.t = State.t + Gen.tstep;

    fprintf('\rIterations = %d | t = %.3f s', it, State.t);
    drawnow;
end

%% - Plotting
Plotter_simulator(Gen,Plotting,Storage,Wells);

