clc
clear all
close all

%% - About
%This reservoir simulator uses a FVM formulation to solve the continuity of
%mass balance equation in 2-D. The boundary conditions are no flow
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

    if State.t >= Gen.t_applied
        Gen.F = Gen.F_2;
    end

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

    State.t = State.t + Gen.tstep;

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
		Storage.Sig_xx(State.step,:) = State.Sig_xx;
        Storage.Sig_yy(State.step,:) = State.Sig_yy;
        Storage.Sig_xy(State.step,:) = State.Sig_xy;
        Storage.e_xx(State.step,:) = State.e_xx;
        Storage.e_yy(State.step,:) = State.e_yy;
        Storage.gamma_xy(State.step,:) = State.gamma_xy;
    end

    fprintf('\rIterations = %d | t = %.3f s', it, State.t);
    drawnow;
end

%% - Plotting
Plotter_simulator(Gen,Plotting,Storage,Wells);

%% - Compare to Mandel
if Gen.Mandel == 1
    %Drained Young's Modulus
    E = Gen.E(1,1);
    %Drained poisson's ratio
    vd = Gen.v(1,1);
    %Size of quarter domain
    a = Gen.Lx;
    b = Gen.Ly;
    %Fluid compressibility
    cf = Flow.cf(1,1);
    %Permeability
    k = Flow.kx0(1,1);
    %Fluid viscosity
    muf = Flow.muf(1,1);
    %Porosity
    phi = Flow.phi0(1,1);
    %Grain bulk modulus
    Ks = Flow.ks;
    %Load applied [Pa]
    F = Gen.F;


    %Dimensionless x for simulator
    x_d = Storage.x/a;
    %Find middle cells
    [y_mid,~] = min(abs(Storage.y - Gen.Ly/2));
    idx = find(Storage.y == y_mid+Gen.Ly/2);

    %Show match to Mandel solution
    fh = figure;
    ax = axes;
    set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
    set(ax,'ActivePositionProperty','position')
    set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
    hold on

    for i = 8:round(length(Storage.TStorage)/7):length(Storage.TStorage)
        %Real time; dimensionless time td = c*t/a^2
        t = Storage.TStorage(i) - Gen.t_applied;
        [Mandel] = Mandel_Comparison(a,b,cf,E,F,k,Ks,muf,phi,t,vd);
        P_d = Mandel.p0;
        h1 = plot(x_d(idx), (Storage.P(i,idx)-Storage.P(1,idx))./P_d,'k','LineWidth',1);
        h2 = plot(Mandel.x, Mandel.Pd,'r--','LineWidth',Plotting.lwidth_1col);
    end

    xlab = xlabel('Distance, $$x/a$$');
    ylab = ylabel('Pressure, $$(P-P(t=0))/P_0$$');
    set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    set(fh, 'Color','white')
    set(gca, 'Box','off', 'TickDir','out');
    lgd = legend([h1 h2],'Simulation','Analytical Soln.','Location','best');
    set(lgd,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    legend box off
    exportgraphics(fh,'Mandel_matlab.pdf','ContentType','vector')


    %Show Mandel-Cryer effect
    fh = figure;
    ax = axes;
    set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
    set(ax,'ActivePositionProperty','position')
    set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');

    h1 = plot((Storage.TStorage - Gen.t_applied)/Mandel.tnorm, Storage.P(:,idx(1))./(P_d + Gen.Pi),'k','LineWidth',1);

    xlab = xlabel('Time, $$t\cdot c_v/a^2$$');
    ylab = ylabel('Pressure, $$P/(P_0+P(t=0))$$');
    set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    set(fh, 'Color','white')
    set(gca, 'Box','off', 'TickDir','out');
    ylim([0 1])

    exportgraphics(fh,'Mandel_pressure_matlab.pdf','ContentType','vector')

    %Show match to Mandel solution, displacement x
    fh = figure;
    ax = axes;
    set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
    set(ax,'ActivePositionProperty','position')
    set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
    hold on

    for i = 8:round(length(Storage.TStorage)/7):length(Storage.TStorage)
        %Real time; dimensionless time td = c*t/a^2
        t = Storage.TStorage(i) - Gen.t_applied;
        [Mandel] = Mandel_Comparison(a,b,cf,E,F,k,Ks,muf,phi,t,vd);
        h1 = plot(x_d(idx), cumsum(Storage.e_xx(i,idx)*Gen.dx),'k','LineWidth',1);
        h2 = plot(Mandel.x, Mandel.ux,'r--','LineWidth',Plotting.lwidth_1col);
    end

    xlab = xlabel('Distance, $$x/a$$');
    ylab = ylabel('x-displacement [m]');
    set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    set(fh, 'Color','white')
    set(gca, 'Box','off', 'TickDir','out');
    lgd = legend([h1 h2],'Simulation','Analytical Soln.','Location','best');
    set(lgd,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    legend box off
    exportgraphics(fh,'Mandel_disp_x_matlab.pdf','ContentType','vector')

    %Show match to Mandel solution, displacement y
    fh = figure;
    ax = axes;
    set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
    set(ax,'ActivePositionProperty','position')
    set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
    hold on

    for i = 8:round(length(Storage.TStorage)/7):length(Storage.TStorage)
        %Real time; dimensionless time td = c*t/a^2
        t = Storage.TStorage(i) - Gen.t_applied;
        [Mandel] = Mandel_Comparison(a,b,cf,E,F,k,Ks,muf,phi,t,vd);
        h1 = plot(x_d(idx), cumsum(Storage.e_yy(i,idx)*Gen.dx),'k','LineWidth',1);
        h2 = plot(Mandel.x, Mandel.uy,'r--','LineWidth',Plotting.lwidth_1col);
    end

    xlab = xlabel('Distance, $$x/a$$');
    ylab = ylabel('y-displacement [m]');
    set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    set(fh, 'Color','white')
    set(gca, 'Box','off', 'TickDir','out');
    lgd = legend([h1 h2],'Simulation','Analytical Soln.','Location','best');
    set(lgd,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    legend box off
    exportgraphics(fh,'Mandel_disp_y_matlab.pdf','ContentType','vector')
end