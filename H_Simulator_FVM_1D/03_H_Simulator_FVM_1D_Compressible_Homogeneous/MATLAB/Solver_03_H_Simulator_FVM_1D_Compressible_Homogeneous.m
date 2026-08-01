clc
clear all
close all

%% - About
%This reservoir simulator uses a FVM formulation to solve the continuity of
%mass balance equation in 1-D. The boundary conditions are fixed pressure
%at the edges (Dirichlet). It uses a constant, homogeneous permeability and porosity.
%There is no gravity and the simulator is single phase. The fluid is
%considered to be slightly compressible. 

%Barnaby Fryer, 2026

%% - Inputs
[Flow, Gen, Plotting, State, Storage] = InputData();

%% - Run Simulation
while State.t < Gen.tf

%% - Solve For Pressure
%Solve for new pressure
[State] = FIMPressure1D_1Phase(Flow,Gen,State);

%Update to new time
State.t = State.t + Gen.tstep;

    %Store results
    if ~isempty(find(State.t == Storage.TStorage,1))
        State.step = State.step + 1;
        Storage.P(State.step,:) = State.P;
        Storage.flux(State.step,:) = State.flux;
    end


end

%% - Validation
%Semi infinite solution, valid at early times
alpha = Flow.kx/(Flow.phi*Flow.muf*Flow.cf);
P_an = 1e5 + (Gen.PL - 1e5)*erfc(Storage.x./(2*sqrt(alpha * Storage.TStorage(2)')));

%Error calculation
[Ep,ind] = max(abs((Storage.P(2,:) - P_an)./(Gen.PL - Gen.PR)));

%% - Plotting pressure
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
% hold on

% xlab = xlabel('Position, $$x$$ [m]');
% ylab = ylabel('Pressure, $$P_{\mathrm{p}}$$ [MPa]');
% set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
% set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
% set(fh, 'Color','white')
% set(gca, 'Box','off', 'TickDir','out');

for i = 1:10:length(Storage.TStorage)

    cla(ax)          % Clear previous lines

    hold(ax,'on')
    plot(ax, Storage.x, Storage.P(i,:)/1e6, ...
        'k-', 'LineWidth', Plotting.lwidth_1col);

    xlab = xlabel('Position, $$x$$ [m]');
    ylab = ylabel('Pressure, $$P_{\mathrm{p}}$$ [MPa]');
    set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    set(fh, 'Color','white')
    set(ax,'Box','off','TickDir','out')
    ylim([0 Gen.PL/1e6])

    % plot(Storage.x,Storage.P(i,:)/1e6, 'k-','LineWidth',Plotting.lwidth_1col);
    P_an = 1e5 + (Gen.PL - 1e5)*erfc(Storage.x./(2*sqrt(alpha * Storage.TStorage(i)')));
    plot(ax,Storage.x,P_an(1,:)/1e6, 'r--','LineWidth',Plotting.lwidth_1col);
    hold(ax,'off')

    lgd = legend('Simulation','Analytical Soln.');
    set(lgd,'Interpreter','latex','fontsize',Plotting.fsize_1col)
    legend box off

    drawnow
    % exportgraphics(fh,['Readme_images/sim03_pressure_matlab_',num2str(i),'.pdf'],'ContentType','vector')
    exportgraphics(ax,...
    sprintf('Readme_images/frame_%04d.png',i),...
    'Resolution',200)
end


%% - Convert to gif

files = dir('Readme_images/frame_*.png');

for k = 1:length(files)
    img = imread(fullfile(files(k).folder,files(k).name));
    [A,map] = rgb2ind(img,256);

    if k == 1
        imwrite(A,map,'Readme_images/sim03_pressure.gif',...
            'gif','LoopCount',Inf,'DelayTime',0.08);
    else
        imwrite(A,map,'Readme_images/sim03_pressure.gif',...
            'gif','WriteMode','append','DelayTime',0.08);
    end
end

%% - Plotting flux
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');
hold on
for i = 2
    plot(Storage.x,Storage.flux(i,:), 'k-','LineWidth',Plotting.lwidth_1col);
end
xlab = xlabel('Position, $$x$$ [m]');
ylab = ylabel('Mass rate [kg/s]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out', 'YScale', 'log');








