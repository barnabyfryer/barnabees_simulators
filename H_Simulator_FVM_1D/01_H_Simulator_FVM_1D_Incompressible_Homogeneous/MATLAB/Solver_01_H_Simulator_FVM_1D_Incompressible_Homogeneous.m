clc
clear all
close all

%% - About
%This reservoir simulator uses a FVM formulation to solve the continuity of
%mass balance equation in 1-D. The boundary conditions are fixed pressure
%at the edges (Dirichlet). It uses a constant, homogeneous permeability.
%There is no gravity and the simulator is single phase. The fluid is
%considered to be incompressible. The formulation is therefore steady
%state.

%Barnaby Fryer, 2026

%% - Inputs
%Number of cells for discretization [-]
Nx = 100;
%Viscosity [Pa sec]
mu = 1;
%Permeability [m^2]
k = 1;
%Length [m]
Lx = 1;
%Fixed boundary pressure on left [Pa]
PL = 10e6;
%Fixed boundary pressure on right [Pa]
PR = 0;

%% - Plotting inputs
Plotting.lwidth_1col = 0.75;
Plotting.Position_1col_matrix = [4.5 3 6 4.5];
Plotting.fsize_1col = 7;

%% - Preparation
%Predefine memory for Jacobian
A = zeros(Nx,Nx);
%Define cell size
dx = Lx/Nx;
%Define grid position
x = linspace(dx/2,Nx*dx-dx/2,Nx);
%Define edge positions
x_edge = linspace(0,Lx,Nx+1);
%Predefine memory for flow through each cell
q = zeros(Nx,1);
%Predefine memory for velocity across cell boundary
ux = zeros(1,Nx+1);

%% - Build transmissivities
%Interior transmissibility
HTx = k/(mu*dx);
%Boundary transmissibility
HTx_boundary = k/(mu*dx/2);

%% - Build Jacobian
for j = 1:Nx
    %Add cell transmissivity
    if j > 1
        A(j,j) = A(j,j) + HTx;
        A(j,j-1) = A(j,j-1) - HTx;
    %Deal with first cell
    else
        A(j,j) = A(j,j) + HTx_boundary;
        q(j,1) = q(j,1) + PL*HTx_boundary;
    end
    %Add cell transmissivity
    if j < Nx
        A(j,j) = A(j,j) + HTx;
        A(j,j+1) = A(j,j+1) - HTx;
    %Deal with last cell
    else
        A(j,j) = A(j,j) + HTx_boundary;
        q(j,1) = q(j,1) + PR*HTx_boundary;
    end
end

%% - Solve for pressure in each cell
P = A\q;

%% - Find velocity between each cell
%Velocity
for j = 2:Nx
    ux(:,j) = (P(j-1)-P(j)).*HTx;
end
%Find velocity at edges
ux(:,1) = (PL - P(1)).*HTx_boundary;
ux(:,Nx+1) = (P(Nx) - PR).*HTx_boundary;

%% - Find error
P_an = (PL + (PR-PL)*x/Lx);
ux_an = -k*(PR-PL)/(mu*Lx)*(ux*0+1);

Ep = max(abs((P' - P_an)./P_an));
Eq = max(abs((ux - ux_an)./ux_an));

%% - Plotting pressure
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
plot(x,P/1e6, 'k-','LineWidth',Plotting.lwidth_1col);
plot(x,P_an/1e6, 'r--','LineWidth',Plotting.lwidth_1col);
xlab = xlabel('$$x$$ [m]');
ylab = ylabel('$$P_{\mathrm{p}}$$ [MPa]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
lgd = legend('Simulation','Analytical Soln.');
set(lgd,'Interpreter','latex','fontsize',Plotting.fsize_1col)
legend box off

%% - Plotting velocity
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
plot(x_edge,ux, 'k-','LineWidth',Plotting.lwidth_1col);
plot(x_edge,ux_an 'r--','LineWidth',Plotting.lwidth_1col);
xlab = xlabel('$$x$$ [m]');
ylab = ylabel('Darcy Flux $$q_{x}$$ [m/s]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
lgd = legend('Simulation','Analytical Soln.');
set(lgd,'Interpreter','latex','fontsize',Plotting.fsize_1col)
legend box off

