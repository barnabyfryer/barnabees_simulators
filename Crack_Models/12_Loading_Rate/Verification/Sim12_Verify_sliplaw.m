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
filename = 'G21_Fig5a_slip.xlsx';

sheets_Fig5a = sheetnames(filename);

Digitized_G21_Fig5a = struct([]);

for k = 1:length(sheets_Fig5a)

    data = readmatrix(filename,'Sheet',sheets_Fig5a{k});

    Digitized_G21_Fig5a(k).name = sheets_Fig5a{k};
    Digitized_G21_Fig5a(k).l_over_lb = data(:,1);
    Digitized_G21_Fig5a(k).vr_over_v0 = data(:,2);

    token = regexp(sheets_Fig5a{k}, 'f0=([-+]?\d*\.?\d+)', 'tokens', 'once');
    Digitized_G21_Fig5a(k).df0_over_b = str2double(token{1});
end


filename = 'G21_Fig5b_slip.xlsx';

sheets_Fig5b = sheetnames(filename);

Digitized_G21_Fig5b = struct([]);

for k = 1:length(sheets_Fig5b)

    data = readmatrix(filename,'Sheet',sheets_Fig5b{k});

    Digitized_G21_Fig5b(k).name = sheets_Fig5b{k};
    Digitized_G21_Fig5b(k).l_over_lb = data(:,1);
    Digitized_G21_Fig5b(k).vr_over_v0 = data(:,2);

    token = regexp(sheets_Fig5b{k}, 'f0=([-+]?\d*\.?\d+)', 'tokens', 'once');
    Digitized_G21_Fig5b(k).df0_over_b = str2double(token{1});
end


filename = 'G21_Fig5c_slip.xlsx';

sheets_Fig5c = sheetnames(filename);

Digitized_G21_Fig5c = struct([]);

for k = 1:length(sheets_Fig5c)

    data = readmatrix(filename,'Sheet',sheets_Fig5c{k});

    Digitized_G21_Fig5c(k).name = sheets_Fig5c{k};
    Digitized_G21_Fig5c(k).l_over_lb = data(:,1);
    Digitized_G21_Fig5c(k).vr_over_v0 = data(:,2);

    token = regexp(sheets_Fig5c{k}, 'f0=([-+]?\d*\.?\d+)', 'tokens', 'once');
    Digitized_G21_Fig5c(k).df0_over_b = str2double(token{1});
end


%% - Read in MATLAB data

folder = 'MATLAB_outputs';

files_MATLAB = dir(fullfile(folder,'*.csv'));

ResultsMATLAB = struct([]);

for k = 1:length(files_MATLAB)

    filename = files_MATLAB(k).name;
    filepath = fullfile(folder,filename);

    %% Read file as text to extract header information

    fid = fopen(filepath,'r');

    headerLines = {};
    lineNumber = 0;

    while true
        line = fgetl(fid);
        lineNumber = lineNumber + 1;

        if ~ischar(line)
            break
        end

        if startsWith(line,'#')
            headerLines{end+1} = line;
        else
            break
        end
    end

    fclose(fid);

    %% Read numerical data
    %
    % lineNumber is now the line containing the column names,
    % so skip it

    data = readmatrix(filepath,...
        'FileType','text',...
        'Delimiter',',',...
        'NumHeaderLines',lineNumber,...
        'ConsecutiveDelimitersRule','join');

    %% Store filename

    ResultsMATLAB(k).filename = filename;

    %% Extract parameters from header

    ResultsMATLAB(k).a_over_b = NaN;
    ResultsMATLAB(k).Df       = NaN;
    ResultsMATLAB(k).rs_type  = '';
    ResultsMATLAB(k).V0       = NaN;
    ResultsMATLAB(k).DT       = NaN;

    for i = 1:length(headerLines)

        line = headerLines{i};

        if contains(line,'a_over_b')
            ResultsMATLAB(k).a_over_b = sscanf(line,'# a_over_b = %f');

        elseif contains(line,'Delta_f0_over_b')
            ResultsMATLAB(k).Df = sscanf(line,'# Delta_f0_over_b = %f');

        elseif contains(line,'rs_type')
            ResultsMATLAB(k).rs_type = strtrim(extractAfter(line,'='));

        elseif contains(line,'V0_over_Vs')
            ResultsMATLAB(k).V0 = sscanf(line,'# V0_over_Vs = %f');

        elseif contains(line,'Delta_T')
            ResultsMATLAB(k).DT = sscanf(line,'# Delta_T = %f');

        elseif contains(line,'df_dts')
            ResultsMATLAB(k).df_dts = sscanf(line,'# df_dts = %f');

        end

    end

    %% Store solution

    ResultsMATLAB(k).t_over_ts  = data(:,1);
    ResultsMATLAB(k).l_over_lb  = data(:,2);
    ResultsMATLAB(k).vr_over_cs = data(:,3);

    % optional fourth column
    if size(data,2) >= 4
        ResultsMATLAB(k).Veff_over_V0 = data(:,4);
        ResultsMATLAB(k).df0_over_b = data(:,5);
    end

end

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
        'Dfi_(.*?)_' ...
        'dDfdt_(.*?)_' ...
        'V0_(.*?)_' ...
        'ab_(.*?)\.csv'],...
        'tokens','once');

    %% Store filename

    ResultsPython(k).filename = filename;

    %% Store parameters

    ResultsPython(k).rs_type = tokens{1};

    ResultsPython(k).Df = decodeNumber(tokens{2});
    ResultsPython(k).dDfdt = decodeNumber(tokens{3});
    ResultsPython(k).V0 = decodeNumber(tokens{4});
    ResultsPython(k).ab = decodeNumber(tokens{5});

    %% Store solution

    ResultsPython(k).t_over_ts  = data(:,1);
    ResultsPython(k).l_over_lb  = data(:,2);
    ResultsPython(k).vr_over_cs = data(:,3);

end

%% - Plot data for crack length, Fig. S5a, a/b = 0.9, DT = 0, slip
fh = figure;

ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');

hold on

for i = 1:length(files_MATLAB)
    if strcmp(ResultsMATLAB(i).rs_type,'slip')
        if ResultsMATLAB(i).a_over_b == 0.9
            plot(ResultsMATLAB(i).l_over_lb,ResultsMATLAB(i).vr_over_cs/ResultsMATLAB(i).V0,'-','LineWidth',Plotting.lwidth_1col,'Color','k')
        end
    end
end

for i = 1:length(files_Python)
    if strcmp(ResultsPython(i).rs_type,'slip')
        if ResultsPython(i).ab == 0.9
            plot(ResultsPython(i).l_over_lb,ResultsPython(i).vr_over_cs/ResultsPython(i).V0,'--','LineWidth',Plotting.lwidth_1col,'Color','m')
        end
    end
end

for k = 1:length(sheets_Fig5a)
    plot(Digitized_G21_Fig5a(k).l_over_lb,Digitized_G21_Fig5a(k).vr_over_v0,'r.')
end

xlab = xlabel('$$\ell/\ell_{b}$$ [-]');
ylab = ylabel('$$v_{\mathrm{r}}/v_0$$ [-]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out', 'XScale','log', 'YScale','log');

axis([1e-1 1e4 5e-1 2e10])

exportgraphics(fh,'Fig5a.pdf','ContentType','vector')


%% - Plot data for crack length, Fig. S5b, a/b = 1.0, DT = 0, slip
fh = figure;

ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');

hold on

for i = 1:length(files_MATLAB)
    if strcmp(ResultsMATLAB(i).rs_type,'slip')
        if ResultsMATLAB(i).a_over_b == 1.0
            plot(ResultsMATLAB(i).l_over_lb,ResultsMATLAB(i).vr_over_cs/ResultsMATLAB(i).V0,'-','LineWidth',Plotting.lwidth_1col,'Color','k')
        end
    end
end

for i = 1:length(files_Python)
    if strcmp(ResultsPython(i).rs_type,'slip')
        if ResultsPython(i).ab == 1.0
            plot(ResultsPython(i).l_over_lb,ResultsPython(i).vr_over_cs/ResultsPython(i).V0,'--','LineWidth',Plotting.lwidth_1col,'Color','m')
        end
    end
end

for k = 1:length(sheets_Fig5b)
    plot(Digitized_G21_Fig5b(k).l_over_lb,Digitized_G21_Fig5b(k).vr_over_v0,'r.')
end

xlab = xlabel('$$\ell/\ell_{b}$$ [-]');
ylab = ylabel('$$v_{\mathrm{r}}/v_0$$ [-]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out', 'XScale','log', 'YScale','log');

axis([1e-1 1e4 5e-1 2e10])

exportgraphics(fh,'Fig5b.pdf','ContentType','vector')

%% - Plot data for crack length, Fig. S5c, a/b = 1.1, DT = 0, slip
fh = figure;

ax = axes;
set(ax,'Units','centimeters','Position',Plotting.Position_1col)
set(ax,'ActivePositionProperty','position')
set(ax,'FontSize',Plotting.fsize_1col,'TickLabelInterpreter','latex');

hold on

for i = 1:length(files_MATLAB)
    if strcmp(ResultsMATLAB(i).rs_type,'slip')
        if ResultsMATLAB(i).a_over_b == 1.1
            plot(ResultsMATLAB(i).l_over_lb,ResultsMATLAB(i).vr_over_cs/ResultsMATLAB(i).V0,'-','LineWidth',Plotting.lwidth_1col,'Color','k')
        end
    end
end

for i = 1:length(files_Python)
    if strcmp(ResultsPython(i).rs_type,'slip')
        if ResultsPython(i).ab == 1.1
            plot(ResultsPython(i).l_over_lb,ResultsPython(i).vr_over_cs/ResultsPython(i).V0,'--','LineWidth',Plotting.lwidth_1col,'Color','m')
        end
    end
end

for k = 1:length(sheets_Fig5c)
    plot(Digitized_G21_Fig5c(k).l_over_lb,Digitized_G21_Fig5c(k).vr_over_v0,'r.')
end

xlab = xlabel('$$\ell/\ell_{b}$$ [-]');
ylab = ylabel('$$v_{\mathrm{r}}/v_0$$ [-]');
set(xlab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(ylab,'Interpreter','latex','fontsize',Plotting.fsize_1col)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out', 'XScale','log', 'YScale','log');

axis([1e-1 1e4 5e-1 2e10])

exportgraphics(fh,'Fig5c.pdf','ContentType','vector')





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