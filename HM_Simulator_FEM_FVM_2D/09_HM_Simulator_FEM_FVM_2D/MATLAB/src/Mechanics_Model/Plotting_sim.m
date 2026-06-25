function [] = Plotting_sim(e_vol,F,Gen,Plotting, Pos,Sig)

%% - Unpack
Sigx = reshape(Sig(:,1),Gen.Nx,Gen.Ny);                     %[Nx,Ny]
Sigy = reshape(Sig(:,2),Gen.Nx,Gen.Ny);                     %[Nx,Ny]
Sigxy = reshape(Sig(:,3),Gen.Nx,Gen.Ny);                    %[Nx,Ny]
x = Pos.x;
y = Pos.y;
x_new = Pos.x_new;
y_new = Pos.y_new;
du = Pos.du;
dv = Pos.dv;
nodesx = Gen.nodesx;
nodesy = Gen.nodesy;

lwidth_1col = Plotting.lwidth_1col;
Position_1col_matrix = Plotting.Position_1col_matrix;
fsize_1col = Plotting.fsize_1col;

%% - Find cell centers
%Make x and y vectors for new cell centers
xc = zeros(Gen.Nx,Gen.Ny);                                  %[Nx,Ny]
yc = zeros(Gen.Nx,Gen.Ny);                                  %[Nx,Ny]

%% - Find forces
%x-direction forces
Fx = F(1:2:end-1);                                          %[Nn,1]
Fx = reshape(Fx,Gen.Nx+1,Gen.Ny+1);                         %[Nx+1,Ny+1]
%y-direction forces
Fy = F(2:2:end);                                            %[Nn,1]
Fy = reshape(Fy,Gen.Nx+1,Gen.Ny+1);                         %[Nx+1,Ny+1]

%Average x and y locations of all 4 nodes of element
for i = 1:Gen.Nx
    xc(i,:) = (Pos.x_new(i,1:Gen.Ny) + Pos.x_new(i,2:Gen.Ny+1)...
        + Pos.x_new(i+1,1:Gen.Ny) + Pos.x_new(i+1,2:Gen.Ny+1))/4;   %[Nx,Ny]
    yc(i,:) = (Pos.y_new(i,1:Gen.Ny) + Pos.y_new(i,2:Gen.Ny+1)...
        + Pos.y_new(i+1,1:Gen.Ny) + Pos.y_new(i+1,2:Gen.Ny+1))/4;   %[Nx,Ny]
end

%% - Plot grid change in position
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',fsize_1col,'TickLabelInterpreter','latex');
hold on
for j = 1:Gen.Nx+1   
    plot(x(j,:),y(j,:),'bl-')
    plot(x(j,:),y(j,:),'bl*')
    plot(x_new(j,:),y_new(j,:),'r-')
    plot(x_new(j,:),y_new(j,:),'r*')
end
for k = 1:Gen.Ny+1
    plot(x(:,k),y(:,k),'bl-')
    plot(x_new(:,k),y_new(:,k),'r-')
end
hold off
xlab = xlabel('x-Location [m]');
ylab = ylabel('y-Location [m]');
set(xlab,'Interpreter','latex','fontsize',fsize_1col)
set(ylab,'Interpreter','latex','fontsize',fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

%% - Plot x displacements
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',fsize_1col,'TickLabelInterpreter','latex');
surf(x,y,du)
xlab = xlabel('x-Location [m]');
ylab = ylabel('y-Location [m]');
zlab = zlabel('x-displacement [m]');
set(xlab,'Interpreter','latex','fontsize',fsize_1col)
set(ylab,'Interpreter','latex','fontsize',fsize_1col)
set(zlab,'Interpreter','latex','fontsize',fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

%% - Plot y displacements
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',fsize_1col,'TickLabelInterpreter','latex');
surf(x,y,dv)
xlab = xlabel('x-Location [m]');
ylab = ylabel('y-Location [m]');
zlab = zlabel('y-displacement [m]');
set(xlab,'Interpreter','latex','fontsize',fsize_1col)
set(ylab,'Interpreter','latex','fontsize',fsize_1col)
set(zlab,'Interpreter','latex','fontsize',fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

%% - Plot x direction forces
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',fsize_1col,'TickLabelInterpreter','latex');
surf(x,y,Fx)
xlab = xlabel('x-Location [m]');
ylab = ylabel('y-Location [m]');
zlab = zlabel('x-direction force [N]');
set(xlab,'Interpreter','latex','fontsize',fsize_1col)
set(ylab,'Interpreter','latex','fontsize',fsize_1col)
set(zlab,'Interpreter','latex','fontsize',fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

%% - Plot y direction forces
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',fsize_1col,'TickLabelInterpreter','latex');
surf(x,y,Fy)
xlab = xlabel('x-Location [m]');
ylab = ylabel('y-Location [m]');
zlab = zlabel('y-direction force [N]');
set(xlab,'Interpreter','latex','fontsize',fsize_1col)
set(ylab,'Interpreter','latex','fontsize',fsize_1col)
set(zlab,'Interpreter','latex','fontsize',fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

%% - Plot x direction stress
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',fsize_1col,'TickLabelInterpreter','latex');
surf(xc,yc,Sigx)
xlab = xlabel('x-Location [m]');
ylab = ylabel('y-Location [m]');
zlab = zlabel('$$\sigma_{xx}$$ [Pa]');
set(xlab,'Interpreter','latex','fontsize',fsize_1col)
set(ylab,'Interpreter','latex','fontsize',fsize_1col)
set(zlab,'Interpreter','latex','fontsize',fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

%% - Plot y direction stress
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',fsize_1col,'TickLabelInterpreter','latex');
surf(xc,yc,Sigy)
xlab = xlabel('x-Location [m]');
ylab = ylabel('y-Location [m]');
zlab = zlabel('$$\sigma_{yy}$$ [Pa]');
set(xlab,'Interpreter','latex','fontsize',fsize_1col)
set(ylab,'Interpreter','latex','fontsize',fsize_1col)
set(zlab,'Interpreter','latex','fontsize',fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

%% - Plot shear stresses
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',fsize_1col,'TickLabelInterpreter','latex');
surf(xc,yc,Sigxy)
xlab = xlabel('x-Location [m]');
ylab = ylabel('y-Location [m]');
zlab = zlabel('$$\sigma_{xy}$$ [Pa]');
set(xlab,'Interpreter','latex','fontsize',fsize_1col)
set(ylab,'Interpreter','latex','fontsize',fsize_1col)
set(zlab,'Interpreter','latex','fontsize',fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

%% - Plot shear stresses
fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',fsize_1col,'TickLabelInterpreter','latex');
surf(xc,yc,e_vol)
xlab = xlabel('x-Location [m]');
ylab = ylabel('y-Location [m]');
zlab = zlabel('$$\epsilon_{v}$$ [-]');
set(xlab,'Interpreter','latex','fontsize',fsize_1col)
set(ylab,'Interpreter','latex','fontsize',fsize_1col)
set(zlab,'Interpreter','latex','fontsize',fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

%% - Plot fixed nodes
nodesfixed = zeros(Gen.Nn,1);
nodesfixed(nodesx) = 1;
nodesfixed(nodesy) = nodesfixed(nodesy) + 2;
nodesfixed = reshape(nodesfixed,Gen.Nx+1,Gen.Ny+1);

fh = figure;
ax = axes;
set(ax,'Units','centimeters','Position',Position_1col_matrix)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',fsize_1col,'TickLabelInterpreter','latex');
surf(x,y,nodesfixed)
xlab = xlabel('x-Location [m]');
ylab = ylabel('y-Location [m]');
zlab = title('Fixed nodes');
set(xlab,'Interpreter','latex','fontsize',fsize_1col)
set(ylab,'Interpreter','latex','fontsize',fsize_1col)
set(zlab,'Interpreter','latex','fontsize',fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');
colorbar('Ticks',[0,1,2,3],...
         'TickLabels',{'None','x-direction','y-direction','Both'})
end