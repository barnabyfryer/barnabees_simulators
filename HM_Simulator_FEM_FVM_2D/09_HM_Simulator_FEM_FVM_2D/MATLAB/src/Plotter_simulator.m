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

%% - Plotting sigma_xx
sig_xx = reshape(Storage.Sig_xx(end,:)+Gen.biot*Storage.P(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,sig_xx/1e6)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = 'Total stress change, $$\Delta S_{xx}$$ [MPa]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

%% - Plotting sigma_yy
sig_yy = reshape(Storage.Sig_yy(end,:)+Gen.biot*Storage.P(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,sig_yy/1e6)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = 'Total stress change, $$\Delta S_{yy}$$ [MPa]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

%% - Plotting sigma_xy
sig_xy = reshape(Storage.Sig_xy(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,sig_xy/1e6)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = '$$\Delta S_{xy}$$ [MPa]';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

%% - Plotting e_xx
e_xx = reshape(Storage.e_xx(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,e_xx)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = '$$e_{xx}$$';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

%% - Plotting e_yy
e_yy = reshape(Storage.e_yy(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,e_yy)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = '$$e_{yy}$$';
set(c.Label, 'Interpreter', 'latex', 'FontSize', Plotting.fsize_1col);
c.FontSize = Plotting.fsize_1col;

%% - Plotting gamma_xy
gamma_xy = reshape(Storage.gamma_xy(end,:),Gen.Nx,Gen.Ny)';
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
imagesc(Storage.x,Storage.y,gamma_xy)
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Position, $$y$$ [m]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
c = colorbar;
c.Label.String = '$$\gamma_{xy}$$';
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

%% - Make gif

% Pressure range
P_all = Storage.P/1e6;
P_lim = [0, max(P_all(:))];

% Stress range
S_all = Storage.Sig_xx + Gen.biot*Storage.P;
S_lim = [min(S_all(:)), max(S_all(:))]/1e6;
S_lim = [-2 10];

fh = figure;
set(fh,'Color','white')

t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');


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


% Right: stress

ax2 = nexttile;

S_plot = reshape(Storage.Sig_xx(1,:) + Gen.biot*Storage.P(1,:),...
    Gen.Nx,Gen.Ny)';

hp = imagesc(ax2,Storage.x,Storage.y,S_plot/1e6);
clim(ax2,S_lim)

axis(ax2,'image')

set(ax2,...
    'Box','off',...
    'TickDir','out',...
    'FontSize',Plotting.fsize_1col,...
    'TickLabelInterpreter','latex')

xlabel(ax2,'Position, $$x$$ [m]','Interpreter','latex')
ylabel(ax2,'Position, $$y$$ [m]','Interpreter','latex')

c2 = colorbar(ax2);
c2.Label.String = 'Total stress change, $$\Delta S_{xx}$$ [MPa]';
set(c2.Label,'Interpreter','latex',...
    'FontSize',Plotting.fsize_1col)


% Animation loop

for i = 1:length(Storage.TStorage)

    % Update pressure
    hpL.CData = reshape(Storage.P(i,:),Gen.Nx,Gen.Ny)'/1e6;

    % Update stress
    hp.CData = reshape(Storage.Sig_xx(i,:) + Gen.biot*Storage.P(i,:),...
        Gen.Nx,Gen.Ny)'/1e6;

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
        imwrite(A,map,'Readme_images/sim09_pressure_sxx.gif',...
            'gif','LoopCount',Inf,'DelayTime',0.08);
    else
        imwrite(A,map,'Readme_images/sim09_pressure_sxx.gif',...
            'gif','WriteMode','append','DelayTime',0.08);
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