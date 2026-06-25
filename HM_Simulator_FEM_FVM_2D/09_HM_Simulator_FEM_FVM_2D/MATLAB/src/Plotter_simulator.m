function [] = Plotter_simulator(Gen,Plotting,Storage,Wells)

%% - Basic calculations
Dist_well = ((Storage.x - Wells.xP).^2 + (Storage.y - Wells.yP).^2).^0.5;
% Colormap
cmap = parula(length(Storage.TStorage));

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
c.Label.String = 'Log10 Permeability (x), $$k_x$$ [m$$^2$$]';
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
c.Label.String = 'Log10 Permeability (y), $$k_y$$ [m$$^2$$]';
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

%% - Plotting volumetric strain
e_vol = reshape(Storage.e_vol(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,e_vol)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = 'Volumetric strain, $$\epsilon_v$$ [-]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

%% - Plotting volumetric strain error
errP = reshape(Storage.errP(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,errP)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = 'Pressure error [-]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

%% - Plotting volumetric strain error
erre = reshape(Storage.erre(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,erre)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = 'Volumetric strain error [-]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

%% - Plotting pressure versus distance for various times
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
for i = 1:length(Storage.TStorage)
    plot(Dist_well,Storage.P(i,:)/1e6,'.','Color',cmap(i,:))
end
xlab = xlabel('Distance [m]');
ylab = ylabel('Pressure, $$P$$ [MPa]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

% Colorbar
colormap(cmap)
cb = colorbar;
cb.Label.String = 'Time [s]';
cb.Label.Interpreter = 'latex';
% Map colors to actual times
clim([Storage.TStorage(1), Storage.TStorage(end)])


%% - Plotting volumetric strain versus distance for various times
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
for i = 1:length(Storage.TStorage)
    plot(Dist_well,Storage.e_vol(i,:),'.','Color',cmap(i,:))
end
xlab = xlabel('Distance [m]');
ylab = ylabel('Volumetric strain, $$\epsilon_{v}$$ [-]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

% Colorbar
colormap(cmap)
cb = colorbar;
cb.Label.String = 'Time [s]';
cb.Label.Interpreter = 'latex';

% Map colors to actual times
clim([Storage.TStorage(1), Storage.TStorage(end)])

%% - Plotting porosity versus distance for various times
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
for i = 1:length(Storage.TStorage)
    plot(Dist_well,Storage.phi(i,:),'.','Color',cmap(i,:))
end
xlab = xlabel('Distance [m]');
ylab = ylabel('Porosity, $$\phi$$ [-]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

% Colorbar
colormap(cmap)
cb = colorbar;
cb.Label.String = 'Time [s]';
cb.Label.Interpreter = 'latex';

% Map colors to actual times
clim([Storage.TStorage(1), Storage.TStorage(end)])

%% - Sanity check orientation of vectors
% x = reshape(Storage.x,Gen.Nx,Gen.Ny)';
% fh = figure;
% ax = axes;
% set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
% set(ax,'ActivePositionProperty','position')
% set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
% imagesc(Storage.x,Storage.y,x)
% xlab = xlabel('Position, $$x$$ [m]');
% ylab = ylabel('Position, $$y$$ [m]');
% set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
% set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
% set(fh, 'Color','white')
% set(gca, 'Box','off', 'TickDir','out');
% c = colorbar;
% c.Label.String = 'x-Position, $$x$$ [m]';
% set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
% c.FontSize = Plotting.fsize_1col;
% 
% y = reshape(Storage.y,Gen.Nx,Gen.Ny)';
% fh = figure;
% ax = axes;
% set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
% set(ax,'ActivePositionProperty','position')
% set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
% imagesc(Storage.x,Storage.y,y)
% xlab = xlabel('Position, $$x$$ [m]');
% ylab = ylabel('Position, $$y$$ [m]');
% set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
% set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
% set(fh, 'Color','white')
% set(gca, 'Box','off', 'TickDir','out');
% c = colorbar;
% c.Label.String = 'y-Position, $$y$$ [m]';
% set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
% c.FontSize = Plotting.fsize_1col;




end