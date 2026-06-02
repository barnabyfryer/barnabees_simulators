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
Nx = 10;
%Viscosity [Pa sec]
mu = 1;
%Permeability, can be heterogeneous but requires redefining k [m^2]
k = 1;
%Length [m]
Lx = 1;
%Fixed boundary pressure on left [Pa]
PL = 10e6;
%Fixed boundary pressure on right [Pa]
PR = 0;

%% - Plotting inputs
Plotting.lwidth_1col = 0.75;
Plotting.Position_1col_matrix = [2.2 1.8 6 4.5];
Plotting.fsize_1col = 7;

%% - Preparation
%Predefine memory for transmissivity at each cell edge
HTx = zeros(1,Nx+1);
%Predefine memory for Jacobian
A = zeros(Nx,Nx);
%Define cell size
dx = Lx/Nx;
%Define grid position
x = linspace(dx/2,Nx*dx-dx/2,Nx);
%Define edge positions
x_edge = linspace(0,Lx,Nx+1);
%Define permeability of each cell
k = ones(1,Nx)*k;
%Define flow through each cell
q = zeros(1*Nx,1);
%Predefine memory for velocity across cell boundary
ux = zeros(1,Nx+1);

%% - Build transmissivities
%Find harmonic average between cells for permeability
kHx = k(1:end-1).*k(2:end)./(k(1:end-1)+k(2:end));
%Calculate transmissivity
HTx(:,2:Nx) = kHx/(mu*dx/2);
%Deal with edge cell boundaries
HTx(:,1) = k(:,1)/(mu*dx/2);
HTx(:,Nx+1) = k(:,Nx)/(mu*dx/2);

%% - Build Jacobian
for j = 1:Nx
    %Add cell transmissivity
    if j > 1
        A(j,j) = A(j,j) + HTx(1,j);
        A(j,j-1) = A(j,j-1) - HTx(1,j);
    %Deal with first cell
    else
        A(j,j) = A(j,j) + HTx(1,j);
        q(j,1) = q(j,1) + PL*HTx(1,j);
    end
    %Add cell transmissivity
    if j < Nx
        A(j,j) = A(j,j) + HTx(1,j+1);
        A(j,j+1) = A(j,j+1) - HTx(1,j+1);
    %Deal with last cell
    else
        A(j,j) = A(j,j) + HTx(1,j+1);
        q(j,1) = q(j,1) + PR*HTx(1,j+1);
    end
end

%% - Solve for pressure in each cell
P = A\q;

%% - Find velocity between each cell
%Velocity
for j = 2:Nx
    ux(:,j) = (P(j-1)-P(j)).*HTx(:,j);
end
%Find velocity at edges
ux(:,1) = (PL - P(1)).*HTx(:,1);
ux(:,Nx+1) = (P(Nx) - PR).*HTx(:,Nx+1);

%% - Plotting pressure
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
plot(x,P/1e6, 'k-','LineWidth',Plotting.lwidth_1col);
xlab = xlabel('$$x$$ [m]');
ylab = ylabel('$$P_{\mathrm{p}}$$ [MPa]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

%% - Plotting velocity
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
plot(x_edge,ux, 'k-','LineWidth',Plotting.lwidth_1col);
xlab = xlabel('$$x$$ [m]');
ylab = ylabel('Darcy Flux $$q_{x}$$ [m/s]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

