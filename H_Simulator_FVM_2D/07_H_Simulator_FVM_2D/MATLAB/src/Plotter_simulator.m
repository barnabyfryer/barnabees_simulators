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

%% - Make gif


% Pressure range
P_all = Storage.P/1e6;
P_lim = [0, max(P_all(:))];

fh = figure;
set(fh,'Color','white')

% Left: pressure

ax1 = nexttile;

P_plot = reshape(Storage.P(1,:),Gen.Nx,Gen.Ny)'/1e6;

hpL = imagesc(ax1,Storage.x,Storage.y,P_plot);
clim(ax1,P_lim)

axis(ax1,'image')

set(ax1,...
    'Box','off',...
    'TickDir','out',...
    'FontSize',Plotting.fsize_1col,...
    'TickLabelInterpreter','latex')

xlabel(ax1,'Position, $$x$$ [m]','Interpreter','latex')
ylabel(ax1,'Position, $$y$$ [m]','Interpreter','latex')

c1 = colorbar(ax1);
c1.Label.String = 'Pressure, $$P_{\mathrm{p}}$$ [MPa]';
set(c1.Label,'Interpreter','latex',...
    'FontSize',Plotting.fsize_1col)

colormap(gray)

% Animation loop

for i = 1:length(Storage.TStorage)

    % Update pressure
    hpL.CData = reshape(Storage.P(i,:)/1e6,Gen.Nx,Gen.Ny)';

    drawnow

    exportgraphics(fh,...
        sprintf('Readme_images/frame_%04d.png',i),...
        'Resolution',200)

end

%Make gif
files = dir('Readme_images/frame_*.png');

for k = 1:length(files)
    img = imread(fullfile(files(k).folder,files(k).name));
    [A,map] = rgb2ind(img,256);

    if k == 1
        imwrite(A,map,'Readme_images/sim07_pressure.gif',...
            'gif','LoopCount',Inf,'DelayTime',0.04);
    else
        imwrite(A,map,'Readme_images/sim07_pressure.gif',...
            'gif','WriteMode','append','DelayTime',0.04);
    end
end

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