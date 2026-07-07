%% - Main file to run
% Script solving the EoM of (Garagash, Phil. Trans. Roy. Soc, 2021)
% Used in Fryer et al., 2026
% Authors: Barnaby Fryer, Mathias Lebihain, Dmitry Garagash, François Passelègue
% Contact: barnaby.fryer@geoazur.unice.fr

clear; clc; close all;

%% -------------------- Plotting --------------------
fsize = 14;
lwidth = 1.4;

%% -------------------- Parameters --------------------
% Friction parameters
Param.a_over_b = 0.55;            %a/b
Param.Delta_f0_over_b = -1;       %Overstress
Param.rs_type = "aging";          %'slip' or 'aging'

% Foreshock and rupture parameters
delta_a_over_L = 1.25e-6 / 0.192e-6;       %Foreshock slip divided by L
Param.V0_over_Vs = 5.1503149729886135e-05; %Ambient sliding velocity
C = 0.3;                                   %Coefficient for hypocentral force

% Approximate minimum of 𝒱 from G21
V_min = 2.138;

% Shear wave speed [m/s], for post processing only
cs = 1800;

%% -------------------- Crack-length domain --------------------
l_ini_over_lb = 3e-5;    % initial crack length
l_fin_over_lb = 1e4;     % final crack length

%% -------------------- Solver tolerances --------------------
% Recommended both 1e-9 for aging law, 1e-10 for slip law
Param.RelTol = 1e-7;
Param.AbsTol = 1e-7;

%% -------------------- Run solver for all foreshocks --------------------
for i = 1:length(delta_a_over_L)
    %% -------------------- Basic calculations --------------------
    % Scaled hypocentral force
    Param.Delta_T = C * delta_a_over_L(i);
    % Scaled ambient rupture velocity
    Param.v0_over_cs = Param.V0_over_Vs;
    % Scaled \bar{v}_0
    Param.bar_v0_over_cs = exp(-Param.Delta_f0_over_b) * Param.v0_over_cs;

    %% -------------------- Solve ODE --------------------
    [t_over_ts, l_over_lb, vr_over_cs, V_eff, reason] = solve_ode_in_l( ...
        l_ini_over_lb, l_fin_over_lb, V_min, Param);


    %% - Export data
    % --- Helper for safe numeric strings
    num2tag = @(x,fmt) strrep(strrep(sprintf(fmt,x),'-','m'),'.','p');

    a_tag  = num2tag(Param.a_over_b,        '%.2f');
    df_tag = num2tag(Param.Delta_f0_over_b, '%.2f');
    DT_tag = num2tag(Param.Delta_T,          '%.3e');
    V0_tag  = num2tag(Param.V0_over_Vs,      '%.3e');

    rs_tag = char(Param.rs_type);  % 'aging' or 'slip'

    outfile = sprintf( ...
        'rupture_ab%s_df%s_%s_DT%s_V0%s.csv', ...
        a_tag, df_tag, rs_tag, DT_tag, V0_tag);



    fid = fopen(strcat('Output/',outfile), 'w');

    % --- Safety check
    if fid == -1
        error('Could not open file for writing.');
    end

    % --- Metadata header
    fprintf(fid, '# Rupture EoM solution (Garagash 2021 / Fryer et al.)\n');
    fprintf(fid, '#\n');
    fprintf(fid, '# a_over_b = %.6g\n', Param.a_over_b);
    fprintf(fid, '# Delta_f0_over_b = %.6g\n', Param.Delta_f0_over_b);
    fprintf(fid, '# rs_type = %s\n', Param.rs_type);
    fprintf(fid, '# V0_over_Vs = %.6e\n', Param.V0_over_Vs);
    fprintf(fid, '# bar_v0_over_cs = %.6e\n', Param.bar_v0_over_cs);
    fprintf(fid, '# Delta_T = %.6e\n', Param.Delta_T);
    fprintf(fid, '# V_min = %.6g\n', V_min);
    fprintf(fid, '# l_ini_over_lb = %.6e\n', l_ini_over_lb);
    fprintf(fid, '# l_fin_over_lb = %.6e\n', l_fin_over_lb);
    fprintf(fid, '# stop_reason = %d\n', reason);
    fprintf(fid, '#\n');

    fprintf(fid, 't_over_ts,l_over_lb,vr_over_cs,V_eff_over_V0\n');

    data = [t_over_ts(:), l_over_lb(:), vr_over_cs(:), V_eff(:)];

    fprintf(fid, '%.8e,%.8e,%.8e,%.8e\n', data.');
    fclose(fid);

    %% -------------------- Display reason --------------------
    switch reason
        case 0
            disp('Integration stopped: rupture velocity dropped below threshold.');
        case 1
            disp('Integration completed: reached final crack length.');
        case 2
            disp('Integration incomplete: solver did not reach final crack length.');
    end


end

%% --------------- Post process to find nucleation length etc ---------------
if reason > 0
    ind_nuc = length(vr_over_cs) - find(flip(vr_over_cs)*cs < 10, 1);
end


%% -------------------- Plot Results --------------------
%Plot crack length in time
fh = figure;
hold on
plot(t_over_ts,l_over_lb,'k','LineWidth',lwidth)
if reason > 0
    plot(t_over_ts(ind_nuc),l_over_lb(ind_nuc),'ko','LineWidth',lwidth,'MarkerFaceColor','w')
end
xlab = xlabel('$$t/t_{\mathrm{s}}$$ [-]');
ylab = ylabel('$$\ell/\ell_b$$ [-]');
set(xlab,'Interpreter','latex','fontsize',fsize)
set(ylab,'Interpreter','latex','fontsize',fsize)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out', 'XScale','log', 'YScale','log');

%Plot rupture velocity in time
fh = figure;
hold on
plot(t_over_ts,vr_over_cs,'k','LineWidth',lwidth)
if reason > 0
    plot(t_over_ts(ind_nuc),vr_over_cs(ind_nuc),'ko','LineWidth',lwidth,'MarkerFaceColor','w')
end
xlab = xlabel('$$t/t_{\mathrm{s}}$$ [-]');
ylab = ylabel('$$v_{\mathrm{r}}/c_{\mathrm{s}}$$ [-]');
set(xlab,'Interpreter','latex','fontsize',fsize)
set(ylab,'Interpreter','latex','fontsize',fsize)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out', 'XScale','log', 'YScale','log');

%Plot slip velocity in time
fh = figure;
hold on
plot(t_over_ts,V_eff,'k','LineWidth',lwidth)
if reason > 0
    plot(t_over_ts(ind_nuc),V_eff(ind_nuc),'ko','LineWidth',lwidth,'MarkerFaceColor','w')
end
xlab = xlabel('$$t/t_{\mathrm{s}}$$ [-]');
ylab = ylabel('$$V_{\mathrm{eff}}/V_0$$ [-]');
set(xlab,'Interpreter','latex','fontsize',fsize)
set(ylab,'Interpreter','latex','fontsize',fsize)
set(fh, 'Color','white')
set(gca, 'Box','off', 'TickDir','out', 'XScale','log', 'YScale','log');
