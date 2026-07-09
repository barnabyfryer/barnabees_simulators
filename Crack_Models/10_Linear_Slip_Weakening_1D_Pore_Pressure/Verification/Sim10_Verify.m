clc
clear all
close all

%% - Plotting inputs
Plotting.Position_3col = [2.2 1.8 5 4.2];
Plotting.Position_1col = [2.2 1.8 7.0 5.5];
Plotting.Position_1col_matrix = [2.2 1.8 6 4.5];
Plotting.Position_1col_cbar = [2.2 1.8 6.0 5.5];
Plotting.Position_1col_inset = [2.2 1.8 2.0 5.5/3];
Plotting.Position_4col_cbar = [2.2 1.8 3.5 3.5];
Plotting.fsize_3col = 7;
Plotting.fsize_1col = 7;
Plotting.fsize_1col_inset = 5;
Plotting.lwidth_1col_inset = 0.5;
Plotting.lwidth_1col = 0.75;
Plotting.lwidth_1col_big = 1;
Plotting.msize_1col = 5;
Plotting.msize_1col_big = 8;
Plotting.cb_width = 0.3;   % cm
Plotting.cb_width2 = 0.2;   % cm
Plotting.gap = 0.2;        % space between axes and colorbar
Plotting.gap2 = -0.6;        % space between axes and colorbar

%% - Read in digitized G&G 12 data
filename = 'G12_cracklength.xlsx';

sheets = sheetnames(filename);

Digitized_length = struct([]);

for k = 1:length(sheets)

    data = readmatrix(filename,'Sheet',sheets{k});

    Digitized_length(k).name = sheets{k};
    Digitized_length(k).time = data(:,1);
    Digitized_length(k).a = data(:,2);

end
 

filename = 'G12_slip.xlsx';

sheets = sheetnames(filename);

Digitized_slip = struct([]);

for k = 1:length(sheets)

    data = readmatrix(filename,'Sheet',sheets{k});

    Digitized_slip(k).name = sheets{k};
    Digitized_slip(k).time = data(:,1);
    Digitized_slip(k).slip = data(:,2);

end

%% - Read in MATLAB data
folder = '..\MATLAB\Output';

files = dir(fullfile(folder,'*.csv'));

Results = struct([]);

for k = 1:length(files)

    filename = files(k).name;
    filepath = fullfile(folder,filename);

    %% Determine file type

    isSSY = contains(filename,'_ssy');

    %% Read data

    data = readmatrix(filepath);

    %% Extract parameters from filename

    tokens = regexp(filename,...
        ['Sim10_dPprod_(.*?)_dtprod_(.*?)_dPinj_(.*?)_' ...
         'taub_(.*?)_fr_(.*?)' ...
         '(_ssy)?\.csv'],...
         'tokens','once');

    dPprod = str2double(strrep(strrep(tokens{1},'m','-'),'p','.'));
    dtprod = str2double(strrep(tokens{2},'p','.'));
    dPinj  = str2double(strrep(tokens{3},'p','.'));
    taub   = str2double(strrep(tokens{4},'p','.'));
    fr     = str2double(strrep(tokens{5},'p','.'));

    %% Store

    Results(k).filename = filename;

    Results(k).dPprod = dPprod;
    Results(k).dtprod = dtprod;
    Results(k).dPinj  = dPinj;
    Results(k).taub   = taub;
    Results(k).fr     = fr;

    if isSSY
        Results(k).time = data(:,1);
        Results(k).a    = data(:,2);
        Results(k).ssy = 1;
    else
        Results(k).time        = data(:,1);
        Results(k).a           = data(:,2);
        Results(k).slip_center = data(:,3);
        Results(k).ssy = 0;
    end

end

%% - Plot data for crack length
fh = figure;

ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');

hold on
for i = 1:length(files)
    if Results(i).ssy == 0
        plot(Results(i).time,Results(i).a,'LineWidth',Plotting.lwidth_1col,'Color','k')
    else
        plot(Results(i).time,Results(i).a,'LineWidth',Plotting.lwidth_1col,'Color',[0.5 0.5 0.5])
    end
end

for i = 1:length(Digitized_length)
    plot(Digitized_length(i).time,Digitized_length(i).a,'r.')
end

xlab = xlabel('Sqrt. of time, $$\tilde{t}^{\frac{1}{2}}$$ [-]');
ylab = ylabel('Crack length, $$\tilde{a}$$ [-]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

axis([0 20 0 7])

%% - Plot data for slip
fh = figure;

ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');

hold on
for i = 1:length(files)
    if Results(i).ssy == 0
        plot(Results(i).time,Results(i).slip_center,'LineWidth',Plotting.lwidth_1col,'Color','k')
    end
end

for i = 1:length(Digitized_slip)
    plot(Digitized_slip(i).time,Digitized_slip(i).slip,'r.')
end


xlab = xlabel('Sqrt. of time, $$\tilde{t}^{\frac{1}{2}}$$ [-]');
ylab = ylabel('Peak slip, $$\tilde{\delta}$$ [-]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out');

axis([0 20 0 2])