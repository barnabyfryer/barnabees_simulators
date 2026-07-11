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
filename = 'G21_FigS5d.xlsx';

sheets_FigS5d = sheetnames(filename);

Digitized_G21_FigS5d = struct([]);

for k = 1:length(sheets_FigS5d)

    data = readmatrix(filename,'Sheet',sheets_FigS5d{k});

    Digitized_G21_FigS5d(k).name = sheets_FigS5d{k};
    Digitized_G21_FigS5d(k).l_over_lb = data(:,1);
    Digitized_G21_FigS5d(k).vr_over_v0 = data(:,2);

    token = regexp(sheets_FigS5d{k}, 'f0=([-+]?\d*\.?\d+)', 'tokens', 'once');
    Digitized_G21_FigS5d(k).df0_over_b = str2double(token{1});
end


%% - Read in MATLAB data
% folder = 'MATLAB_noprod';
% 
% files = dir(fullfile(folder,'*.csv'));
% 
% Results = struct([]);
% 
% for k = 1:length(files)
% 
%     filename = files(k).name;
%     filepath = fullfile(folder,filename);
% 
%     %% Determine file type
% 
%     isSSY = contains(filename,'_ssy');
% 
%     %% Read data
% 
%     data = readmatrix(filepath);
% 
%     %% Extract parameters from filename
% 
%     tokens = regexp(filename,...
%         ['Sim10_dPprod_(.*?)_dtprod_(.*?)_dPinj_(.*?)_' ...
%          'taub_(.*?)_fr_(.*?)' ...
%          '(_ssy)?\.csv'],...
%          'tokens','once');
% 
%     dPprod = str2double(strrep(strrep(tokens{1},'m','-'),'p','.'));
%     dtprod = str2double(strrep(tokens{2},'p','.'));
%     dPinj  = str2double(strrep(tokens{3},'p','.'));
%     taub   = str2double(strrep(tokens{4},'p','.'));
%     fr     = str2double(strrep(tokens{5},'p','.'));
% 
%     %% Store
% 
%     Results(k).filename = filename;
% 
%     Results(k).dPprod = dPprod;
%     Results(k).dtprod = dtprod;
%     Results(k).dPinj  = dPinj;
%     Results(k).taub   = taub;
%     Results(k).fr     = fr;
% 
%     if isSSY
%         Results(k).time = data(:,1);
%         Results(k).a    = data(:,2);
%         Results(k).ssy = 1;
%     else
%         Results(k).time        = data(:,1);
%         Results(k).a           = data(:,2);
%         Results(k).slip_center = data(:,3);
%         Results(k).ssy = 0;
%     end
% 
% end

%% - Read in Python data

folder = 'Python_outputs';

files_Python = dir(fullfile(folder,'*.csv'));

ResultsPython = struct([]);

for k = 1:length(files_Python)

    filename = files_Python(k).name;
    filepath = fullfile(folder,filename);

    %% Read data

    data = readmatrix(filepath,...
        'Delimiter',';',...
        'NumHeaderLines',1);

    %% Extract parameters from filename

    tokens = regexp(filename,...
        ['EoM_' ...
         '(.*?)_' ...
         'Df_(.*?)_' ...
         'DT_(.*?)_' ...
         'V0_(.*?)_' ...
         'ab_(.*?)\.csv'],...
         'tokens','once');

    %% Store filename

    ResultsPython(k).filename = filename;

    %% Store parameters

    ResultsPython(k).rs_type = tokens{1};

    ResultsPython(k).Df = decodeNumber(tokens{2});
    ResultsPython(k).DT = decodeNumber(tokens{3});
    ResultsPython(k).V0 = decodeNumber(tokens{4});
    ResultsPython(k).ab = decodeNumber(tokens{5});

    %% Store solution

    ResultsPython(k).t_over_ts  = data(:,1);
    ResultsPython(k).l_over_lb  = data(:,2);
    ResultsPython(k).vr_over_cs = data(:,3);

end


%% - Plot data for crack length, Fig. S5c, a/b = 0.9, DT = 2, aging
fh = figure;

ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');

hold on

for i = 1:length(files_Python)
    if ResultsPython(i).ab == 0.9
        if ResultsPython(i).DT == 2
            plot(ResultsPython(i).l_over_lb,ResultsPython(i).vr_over_cs/ResultsPython(i).V0,'--','LineWidth',Plotting.lwidth_1col,'Color','m')
        end
    end
end

for k = 1:length(sheets_FigS5d)
    plot(Digitized_G21_FigS5d(k).l_over_lb,Digitized_G21_FigS5d(k).vr_over_v0,'r.')
end

xlab = xlabel('$$\ell/\ell_{b}$$ [-]');
ylab = ylabel('$$v_{\mathrm{r}}/v_0$$ [-]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out', 'XScale','log', 'YScale','log');

axis([1e-1 1e4 5e-1 2e10])




%% - For reading in file data
function x = decodeNumber(str)

% Undo filename formatting:
%   p  -> .
%   em -> e-
%   e  -> e
%   m  -> -

str = strrep(str,'em','e-');
str = strrep(str,'p','.');
str = strrep(str,'m','-');

x = str2double(str);

end