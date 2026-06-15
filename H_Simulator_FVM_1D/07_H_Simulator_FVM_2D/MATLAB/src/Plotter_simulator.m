function [] = Plotter_simulator(Gen,Plotting,Storage)


%% - Plotting pressure
P_plot = reshape(Storage.P(end,:),Gen.Nx,Gen.Ny)'/1e6;
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,P_plot);
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = 'Pressure, $$P_{\mathrm{p}}$$ [MPa]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

%% - Plotting permeability
Kx = reshape(Storage.kx(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,log10(Kx))
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = 'Permeability (x), $$k_x$$ [m$$^2$$]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

Ky = reshape(Storage.ky(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,log10(Ky))
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = 'Permeability (y), $$k_y$$ [m$$^2$$]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

%% - Plotting porosity
phi = reshape(Storage.phi(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,phi)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = 'Porosity, $$\phi$$ [-]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

%% - Plotting flux
flux = reshape(Storage.flux(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,flux)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = 'Mass flux [kg/s]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;




x = reshape(Storage.x,Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,x)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = 'x-Position, $$x$$ [m]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

y = reshape(Storage.y,Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,y)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = 'y-Position, $$y$$ [m]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;




end