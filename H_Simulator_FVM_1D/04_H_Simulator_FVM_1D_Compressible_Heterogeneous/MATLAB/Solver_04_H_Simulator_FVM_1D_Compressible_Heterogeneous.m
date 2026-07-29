clc
clear all
close all

%% - About
%This reservoir simulator uses a FVM formulation to solve the continuity of
%mass balance equation in 1-D. The boundary conditions are fixed pressure
%at the edges (Dirichlet). It uses a constant, heterogeneous permeability and porosity.
%There is no gravity and the simulator is single phase. The fluid is
%considered to be slightly compressible. 

%Barnaby Fryer, 2026

%% - Inputs
[Flow, Gen, Plotting, State, Storage] = InputData();

%% - Run Simulation
while State.t < Gen.tf

%% - Solve For Pressure
[State] = FIMPressure1D_1Phase(Flow,Gen,State);

    %Store results
    if ~isempty(find(State.t == Storage.TStorage,1))
        State.step = State.step + 1;
        Storage.P(State.step,:) = State.P;
        Storage.flux(State.step,:) = State.flux;
    end

State.t = State.t + Gen.tstep;
end

%% - Plotting pressure
%Steady state analytcal solution for two permeability zones
x = linspace(0,Gen.Lx,100);
p_an = Gen.PL - (Gen.PL - Gen.PR)*x/Flow.kx(1)/((Flow.L_k/Flow.kx(1)) + (Gen.Lx - Flow.L_k)/Flow.kx(end));
p_an(x > Flow.L_k) = Gen.PL - (Gen.PL - Gen.PR)*(Flow.L_k/Flow.kx(1) + (x(x > Flow.L_k)-Flow.L_k)/Flow.kx(end))/(Flow.L_k/Flow.kx(1) + (Gen.Lx - Flow.L_k)/Flow.kx(end));


fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
for i = length(Storage.TStorage)
    plot(Storage.x,Storage.P(i,:)/1e6, 'k-','LineWidth',Plotting.lwidth_1col);
end
plot(x,p_an/1e6, 'r--','LineWidth',Plotting.lwidth_1col);
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Pressure, $$P_{\mathrm{p}}$$ [MPa]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
lgd = legend('Simulation','Analytical Soln.');
set(lgd,'Interpreter','latex','fontsize',Plotting.fsize_1col)
legend box off

%% - Plotting permeability
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
plot(Storage.x,Flow.kx, 'k-','LineWidth',Plotting.lwidth_1col);
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Permeability, $$k$$ [m$$^2$$]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out', 'YScale', 'log');

%% - Plotting porosity
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
plot(Storage.x,Flow.phi, 'k-','LineWidth',Plotting.lwidth_1col);
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Porosity, $$\phi$$ [-]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

%% - Plotting flux
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
for i = 1:length(Storage.TStorage)
    plot(Storage.x,Storage.flux(i,:), 'k-','LineWidth',Plotting.lwidth_1col);
end
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Mass flux [kg/s]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out', 'YScale', 'log');








