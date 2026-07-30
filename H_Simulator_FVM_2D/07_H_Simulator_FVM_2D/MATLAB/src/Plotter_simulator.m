function [] = Plotter_simulator(Flow,Gen,Plotting,Storage,Wells)


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

%% - Theis solution
dist_max = ((0 - Wells.xQ(1)).^2 + (0 - Wells.yQ(1)).^2).^(0.5);
dist_dummy = linspace(0,dist_max,100);
t_eval = Storage.TStorage(end);
ct = Flow.cf;
alpha = Flow.kx0(1)/(Flow.phi0(1)*Flow.muf*ct);
P_Theis = Storage.P(1,1) + Flow.muf*Wells.Q(1)/(4*pi*Flow.kx0(1)*Gen.Lz*Flow.Rho0)*expint(dist_dummy.^2/(4*alpha*t_eval));

dist = ((Storage.x - Wells.xQ(1)).^2 + (Storage.y - Wells.yQ(1)).^2).^(0.5);
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
plot(dist,Storage.P(end,:)/1e6,'k.');
plot(dist_dummy,P_Theis/1e6,'r--','LineWidth',Plotting.lwidth_1col);
xlab = xlabel('Radial distanec [m]');
ylab = ylabel('Pressure, $$P$$ [MPa]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
lgd = legend('Simulation','Analytical Soln.','Location','best');
set(lgd,'Interpreter','latex','fontsize',Plotting.fsize_1col)
legend box off

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