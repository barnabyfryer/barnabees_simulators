function [] = Plotter_simulator(Flow,Gen,Plotting,Storage,Wells)


%% - Plotting pressure
% Here we find an "analytical" solution by reaching steady state with one
% fixed pressure well and one fixed injection mass rate. Then, knowing the
% mass rate everywhere, we find the pressure at the next cell. This results
% in an update to permeability, which in turn updates pressure again. We
% keep iterating until convergence on that cell. Then we move to the next
% cell.
% Q = -Wells.Q(1);
% tol = 1e-6;
% P = Wells.P(1);
% for i = 1:length(Storage.x)
%     err = 1;
% 
%     P(i+1) = P(i);      % initial guess
% 
%     while err > tol
% 
%         P_old = P(i+1);
%         State.P = P_old;
%         % permeability evaluated from current pressure guess
%         State = Perm(Flow,State);
%         k = State.kx(1);
%         % Darcy equation
%         P(i+1) = P(i) - Q*Flow.muf*Gen.dx/(Flow.Rho0*k)*(Gen.Ly*Gen.Lz);
% 
%         err = abs(P(i+1)-P_old);
% 
%     end
% end


fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
for i = length(Storage.TStorage)
    plot(Storage.x,Storage.P(i,:)/1e6, 'k-','LineWidth',Plotting.lwidth_1col);
end
% plot(Storage.x,P(1:end-1)/1e6, 'r--','LineWidth',Plotting.lwidth_1col);
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Pressure, $$P_{\mathrm{p}}$$ [MPa]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
lgd = legend('Simulation','Analytical Soln.','Location','best');
set(lgd,'Interpreter','latex','fontsize',Plotting.fsize_1col)
legend box off

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
for i = length(Storage.TStorage)
    plot(Storage.x,Storage.flux(i,:), 'k-','LineWidth',Plotting.lwidth_1col);
end
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Mass rate [kg/s]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out', 'YScale', 'log');



% %% - Plot total mass for validation
% 
% fh = figure;
% ax = axes;
% set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
% set(ax,'ActivePositionProperty','position')
% set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
% hold on
% t = linspace(0,max(Storage.TStorage),100);
% plot(t,t*Wells.Q, 'r--','LineWidth',Plotting.lwidth_1col);
% for i = 1:length(Storage.TStorage)
%     rho = Density(Flow,Storage.P(i,:));
%     State.P = Storage.P(i,:)';
%     phi = PhiCalc(Flow,State);
%     V = Gen.dx*Gen.Ly*Gen.Lz;
%     mass(i) = sum(rho'.*phi.*V);
%     plot(Storage.TStorage(i),mass(i) - mass(1), 'ko','LineWidth',Plotting.lwidth_1col,'MarkerFaceColor','w');
% end
% xlab = xlabel('Time, $$t$$ [s]');
% ylab = ylabel('Mass change [kg]');
% set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
% set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
% set(fh, 'Color','white')
% set(gca, 'Box','off', 'TickDir','out');
% lgd = legend('Analytical Soln.','Simulation','Location','best');
% set(lgd,'Interpreter','latex','fontsize',Plotting.fsize_1col)
% legend box off
% 


%% - Make gif
fh = figure;
set(fh,'Color','white')

t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Left: pressure
ax1 = nexttile;
hpL = plot(ax1,Storage.x,Storage.P(1,:)/1e6,...
    'k-','LineWidth',Plotting.lwidth_1col);
set(ax1,...
    'Box','off',...
    'TickDir','out',...
    'FontSize',Plotting.fsize_1col,...
    'TickLabelInterpreter','latex')
xlabel(ax1,'Position, $$x$$ [m]','Interpreter','latex')
ylabel(ax1,'Pressure, $$P_{\mathrm{p}}$$ [MPa]','Interpreter','latex')
ylim(ax1,[0 Wells.P(1)/1e6])


% Right: permeability
ax2 = nexttile;
hp = plot(ax2,Storage.x,Storage.k(1,:),'k-','LineWidth',Plotting.lwidth_1col);

set(ax2,'Box','off',...
    'TickDir','out',...
    'FontSize',Plotting.fsize_1col,...
    'TickLabelInterpreter','latex')

xlabel(ax2,'Position, $$x$$ [m]','Interpreter','latex')
ylabel(ax2,'Permeability, $$k$$ [m$^2$]','Interpreter','latex')
ylim([1 1.7]*1e-12)

for i = 1:length(Storage.TStorage)

    hpL.YData = Storage.P(i,:)/1e6;

    hp.YData = Storage.k(i,:);

    drawnow

    exportgraphics(fh,...
        sprintf('Readme_images/frame_%04d.png',i),...
        'Resolution',200)

end

files = dir('Readme_images/frame_*.png');

for k = 1:length(files)
    img = imread(fullfile(files(k).folder,files(k).name));
    [A,map] = rgb2ind(img,256);

    if k == 1
        imwrite(A,map,'Readme_images/sim06_pressure_k.gif',...
            'gif','LoopCount',Inf,'DelayTime',0.08);
    else
        imwrite(A,map,'Readme_images/sim06_pressure_k.gif',...
            'gif','WriteMode','append','DelayTime',0.08);
    end
end



end