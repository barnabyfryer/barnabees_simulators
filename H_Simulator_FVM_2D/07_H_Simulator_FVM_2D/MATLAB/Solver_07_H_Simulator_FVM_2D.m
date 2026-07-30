clc
clear all
close all

%% - About
%This reservoir simulator uses a FVM formulation to solve the continuity of
%mass balance equation in 2-D. The boundary conditions are no flow
%at the edges (Neumann). It uses a pressure-dependent, heterogeneous permeability and porosity.
%There is no gravity and the simulator is single phase. The fluid is
%considered to be slightly compressible. 

%Barnaby Fryer, 2026

%% - Inputs
[Flow, Gen, Plotting, State, Storage, Wells] = InputData();

%% - Run Simulation
while State.t < Gen.tf

    %% - Solve For Pressure
    [State] = FIMPressure2D_1Phase(Flow,Gen,State,Wells);

    State.t = State.t + Gen.tstep;

    %Store results
    if ~isempty(find(State.t == Storage.TStorage,1))
        State.step = State.step + 1;
        Storage.P(State.step,:) = State.P;
        Storage.flux(State.step,:) = State.flux;
        Storage.kx(State.step,:) = State.kx;
        Storage.ky(State.step,:) = State.ky;
        Storage.phi(State.step,:) = State.phi;
    end

    
end

%% - Plotting
Plotter_simulator(Flow,Gen,Plotting,Storage,Wells);

