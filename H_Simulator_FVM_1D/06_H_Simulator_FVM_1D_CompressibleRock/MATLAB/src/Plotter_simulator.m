function [] = Plotter_simulator(Plotting,Storage)


%% - Plotting pressure
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
for i = 1:length(Storage.TStorage)
    plot(Storage.x,Storage.P(i,:)/1e6, 'k-','LineWidth',Plotting.lwidth_1col);
end
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Pressure, $$P_{\mathrm{p}}$$ [MPa]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

%% - Plotting permeability
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
plot(Storage.x,Storage.k(end,:), 'k-','LineWidth',Plotting.lwidth_1col);
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
plot(Storage.x,Storage.phi(end,:), 'k-','LineWidth',Plotting.lwidth_1col);
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









end