%%% Figures for, "Increasing Sea Surface Warming Suppresses Primary Sea
%%% Spray Aerosol Number Production" – egusphere-2026-2142

% Justin Hamlin
% 12 August 2026

clear, clc, close all
cd('C:/Users/justi/OneDrive/Documents/GitHub/egusphere-2026-2142/')
addpath('C:/Users/justi/OneDrive/Documents/GitHub/egusphere-2026-2142/')
data = load('soars_polar_final.mat').all_stats;
smops = load('soars_polar_smops_final.mat').all_stats;

bubbles_data = load('soars_polar_bubbles_final.mat').all_stats;

%% Load Colormaps - From RJLIII
cd C:/Users/justi/OneDrive/Documents/0_UCSD/Research/MATLAB/Colormap/slanCM/
load('slanCM_Data.mat')
ColorMaps = slandarerCM;
clear('author','slandarerCM')
CM = flipud(ColorMaps(6).Colors{11});
for i = 1:5
    Colors(i,:) = CM(i*floor(256/5), 1:3);
end
Colors = flipud(Colors);

%%% Figure Standardization
fs = 14; % set font size
lw = 1.5; % set line width
ms = 4; % set marker size
fColor = "k";


%%% ------------------------------------------------------------------ %%%
% pnsd and modal number concentration

figure(1), clf
t = tiledlayout('flow');

nexttile(1, [1 3])
hold on
for i = 1:size(data, 2)
    patch([data{i}.D, fliplr(data{i}.D)], ...
        [data{i}.dN.mean + data{i}.dN.std, ...
        fliplr(data{i}.dN.mean - data{i}.dN.std)], ...
        Colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none')

    plot(data{i}.D, data{i}.dN.mean, Color=Colors(i,:), LineWidth=lw+0.75)
end

nexttile(2)
% combine data into table
yVar = "N_aitken";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.(yVar).Lin_NRMSE, ...
    reg_metrics.(yVar).Lin_p, ...
    reg_metrics.(yVar).Lin_VarExplained, ...
    reg_metrics.(yVar).Quad_NRMSE, ...
    reg_metrics.(yVar).Quad_p, ...
    reg_metrics.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

ylim([0 2500])
ylabel('Aitken N [cm^{-3}]')

nexttile(3)
% combine data into table
yVar = "N_accumulation";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.(yVar).Lin_NRMSE, ...
    reg_metrics.(yVar).Lin_p, ...
    reg_metrics.(yVar).Lin_VarExplained, ...
    reg_metrics.(yVar).Quad_NRMSE, ...
    reg_metrics.(yVar).Quad_p, ...
    reg_metrics.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end
ylim([0 4000])
ylabel('Accumulation N [cm^{-3}]')

nexttile(4)
% combine data into table
yVar = "N_supermicron";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.(yVar).Lin_NRMSE, ...
    reg_metrics.(yVar).Lin_p, ...
    reg_metrics.(yVar).Lin_VarExplained, ...
    reg_metrics.(yVar).Quad_NRMSE, ...
    reg_metrics.(yVar).Quad_p, ...
    reg_metrics.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

ylim([0 8])
ylabel('Supermicron N [cm^{-3}]')

clim([0 25])
colormap(Colors)
c = colorbar;
c.Limits = [0 25];
c.Ticks = [0 5 10 15 20 25];
c.TickLabels = ["0" "5" "10" "15" "20" "25"];
c.Label.String = 'Sea Surface Temperature (SST) [°C]';
c.Label.FontSize = fs+2;
c.TickDirection = 'out';
c.Layout.Tile = 'east';

labels = {'(a)','(b)','(c)','(d)'};
for i = 1:4
    nexttile(i)

    ax = gca;
    ax.FontSize = fs;
    ax.TickDir = 'out';
    if i == 1
        ax.XScale = "log";
        ax.XLim = [10^-2 10^1];
        ax.XLabel.String = 'Dry, Physical Diameter (D_{p, dry}) [\mum]';
        ax.YLabel.String = "PNSD (dN/dlogD_p) [cm^{-3}]";
        ax.YLim = [0 8000];
        ax.TickLength = [0.01 0.025];
        text(0.01, 0.96, labels{i}, 'Units', 'normalized', 'FontSize', fs, 'FontWeight', 'bold')
        
    else
        ax.XScale = "linear";
        text(0.05, 0.96, labels{i}, 'Units', 'normalized', 'FontSize', fs, 'FontWeight', 'bold')
        limits = [0 25];
        ticks = 0:5:25;
        ax.XLim = limits;
        ax.XTick = ticks;
        ax.XTickLabel = string(ticks);
        ax.TickLength = [0.02 0.025];
    end
end

t.XLabel.String = 'Sea Surface Temperature (SST) [°C]';
t.XLabel.FontSize = fs+2;

set(gcf,'Position',[50 50 1150 750],'Color','w') % set standard figure size


%%% ------------------------------------------------------------------ %%%
% Figure: Key Parameters of the PNSD
figure(2), clf
t = tiledlayout('flow');

nexttile(1)
hold on
% combine data into table
yVar = "N";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.(yVar).Lin_NRMSE, ...
    reg_metrics.(yVar).Lin_p, ...
    reg_metrics.(yVar).Lin_VarExplained, ...
    reg_metrics.(yVar).Quad_NRMSE, ...
    reg_metrics.(yVar).Quad_p, ...
    reg_metrics.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

ylim([0 6000])
ylabel('N [cm^{-3}]', 'FontSize', fs)

nexttile(2)
hold on
% combine data into table
yVar = "Dg";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.(yVar).Lin_NRMSE, ...
    reg_metrics.(yVar).Lin_p, ...
    reg_metrics.(yVar).Lin_VarExplained, ...
    reg_metrics.(yVar).Quad_NRMSE, ...
    reg_metrics.(yVar).Quad_p, ...
    reg_metrics.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end
ylim([0.08 0.14])
ylabel('D_g [\mum]', 'FontSize', fs)

nexttile(3)
hold on
% combine data into table
yVar = "sigma_g";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.(yVar).Lin_NRMSE, ...
    reg_metrics.(yVar).Lin_p, ...
    reg_metrics.(yVar).Lin_VarExplained, ...
    reg_metrics.(yVar).Quad_NRMSE, ...
    reg_metrics.(yVar).Quad_p, ...
    reg_metrics.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

ylim([1.8 2.2])
ylabel('\sigma_g', 'FontSize', fs)

nexttile(4)
hold on
% combine data into table
yVar = "D_mode";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.(yVar).Lin_NRMSE, ...
    reg_metrics.(yVar).Lin_p, ...
    reg_metrics.(yVar).Lin_VarExplained, ...
    reg_metrics.(yVar).Quad_NRMSE, ...
    reg_metrics.(yVar).Quad_p, ...
    reg_metrics.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

ylim([0.05 0.2])
ylabel('D_{mode} [\mum]', 'FontSize', fs)

labels = {'(a)','(b)','(c)','(d)'};
for i = 1:4
    nexttile(i)
    text(0.02, 0.96, labels{i}, 'Units', 'normalized', 'FontSize', fs, 'FontWeight', 'bold')

    ax = gca;
    ax.TickLength = [0.02 0.025];
    ax.FontSize = fs;
    ax.TickDir = 'out';
    limits = [0 25];
    ticks = 0:5:25;
    ax.XLim = limits;
    ax.XTick = ticks;
    ax.XTickLabel = string(ticks);

end

t.XLabel.String = 'Sea Surface Temperature (SST) [°C]';
t.XLabel.FontSize = fs;
set(gcf,'Position',[50 50 800 600],'Color','w') % set standard figure size

%%% ------------------------------------------------------------------ %%%
% 6 panel N, S, V and distributions
figure(3), clf
t = tiledlayout(3,2);

nexttile(1)
hold on
% combine data into table
yVar = "N";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

nexttile(2)
hold on
for i = 1:size(data, 2)
    plot(data{i}.D, data{i}.dN.mean, Color=Colors(i,:), LineWidth=lw)
end

nexttile(3)
hold on
% combine data into table
yVar = "S";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.(yVar).Lin_NRMSE, ...
    reg_metrics.(yVar).Lin_p, ...
    reg_metrics.(yVar).Lin_VarExplained, ...
    reg_metrics.(yVar).Quad_NRMSE, ...
    reg_metrics.(yVar).Quad_p, ...
    reg_metrics.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end
nexttile(4)
hold on
for i = 1:size(data, 2)
    plot(data{i}.D, data{i}.dS.mean, Color=Colors(i,:), LineWidth=lw)
end

nexttile(5)
hold on
% combine data into table
yVar = "V";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.(yVar).Lin_NRMSE, ...
    reg_metrics.(yVar).Lin_p, ...
    reg_metrics.(yVar).Lin_VarExplained, ...
    reg_metrics.(yVar).Quad_NRMSE, ...
    reg_metrics.(yVar).Quad_p, ...
    reg_metrics.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

nexttile(6)
hold on
for i = 1:size(data, 2)
    plot(data{i}.D, data{i}.dV.mean, Color=Colors(i,:), LineWidth=lw)
end

labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};
for i = 1:6
    nexttile(i)
    ax = gca;
    if i == 1
        ax.YLabel.String = 'N [cm^{-3}]';
        ax.YLim = [0 6000];
    elseif i == 2
        ax.YLabel.String = 'dN/dlogD_{p, dry} [cm^{-3}]';
        ax.YLim = [0 8000];
    elseif i == 3
        ax.YLabel.String = 'S [\mum^2 cm^{-3}]';
        ax.YLim = [0 800];
    elseif i == 4
        ax.YLabel.String = 'dS/dlogD_{p, dry} [\mum^2 cm^{-3}]';
        ax.YTick = ax.YLim(1):200:ax.YLim(2);
    elseif i == 5
        ax.YLabel.String = 'V [\mum^3 cm^{-3}]';
        ax.YLim = [0 50];
    elseif i == 6
        ax.YLabel.String = 'dV/dlogD_{p, dry} [\mum^3 cm^{-3}]';
    end

    if i == 2 || i == 4 || i == 6
        ax.XScale = 'log';
        ax.XLim = [10^-2 10^1];
        ax.YScale = 'linear';
        ax.YLabel.Position(1) = 0.0028;
    end

    if i == 5
        ax.XLabel.String = 'Sea Surface Temperature (SST) [°C]';
    elseif i == 6
        ax.XLabel.String = 'Dry, Physical Diameter (D_{p, dry}) [\mum]';
    end
    ax.FontSize = fs;
    ax.YLabel.FontSize = fs-1;
    ax.XLabel.FontSize = fs-1;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.025];
    ax.Box = 'off';
    text(0.02, 0.96, labels{i}, 'Units', 'normalized', 'FontSize', fs, 'FontWeight', 'bold')
end

clim([0 25])
colormap(Colors)
c = colorbar;
c.Limits = [0 25];
c.Ticks = [0 5 10 15 20 25];
c.TickLabels = ["0" "5" "10" "15" "20" "25"];
c.Label.String = 'Sea Surface Temperature (SST) [°C]';
c.Label.FontSize = fs;
c.TickDirection = 'out';
c.Layout.Tile = 'east';

set(gcf,'Position',[50 50 1000 800],'Color','w') % set standard figure size

%%% ------------------------------------------------------------------ %%%
% SSA N and subsurface n comparison

figure(4), clf
t = tiledlayout(1,3);

nexttile(1)
hold on
% combine data into table
xVar = "n";
yVar = "N_aitken";
for j = 1:size(data,2)
    tbl(j,:) = table(bubbles_data{j}.(xVar).mean, bubbles_data{j}.(xVar).mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[ssb_reg_metrics.(yVar).Lin_NRMSE, ...
    ssb_reg_metrics.(yVar).Lin_p, ...
    ssb_reg_metrics.(yVar).Lin_VarExplained, ...
    ssb_reg_metrics.(yVar).Quad_NRMSE, ...
    ssb_reg_metrics.(yVar).Quad_p, ...
    ssb_reg_metrics.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
ssb_reg_metrics.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
ssb_reg_metrics.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southeast';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end
ylim([0 2500])

nexttile(2)
hold on
% combine data into table
xVar = "n";
yVar = "N_accumulation";
for j = 1:size(data,2)
    tbl(j,:) = table(bubbles_data{j}.(xVar).mean, bubbles_data{j}.(xVar).mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[ssb_reg_metrics.(yVar).Lin_NRMSE, ...
    ssb_reg_metrics.(yVar).Lin_p, ...
    ssb_reg_metrics.(yVar).Lin_VarExplained, ...
    ssb_reg_metrics.(yVar).Quad_NRMSE, ...
    ssb_reg_metrics.(yVar).Quad_p, ...
    ssb_reg_metrics.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
ssb_reg_metrics.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
ssb_reg_metrics.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;


% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southeast';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end
ylim([0 4000])

nexttile(3)
% combine data into table
xVar = "n";
yVar = "N_supermicron";
for j = 1:size(data,2)
    tbl(j,:) = table(bubbles_data{j}.(xVar).mean, bubbles_data{j}.(xVar).mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[ssb_reg_metrics.(yVar).Lin_NRMSE, ...
    ssb_reg_metrics.(yVar).Lin_p, ...
    ssb_reg_metrics.(yVar).Lin_VarExplained, ...
    ssb_reg_metrics.(yVar).Quad_NRMSE, ...
    ssb_reg_metrics.(yVar).Quad_p, ...
    ssb_reg_metrics.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
ssb_reg_metrics.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
ssb_reg_metrics.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southeast';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end
ylim([0 8])

labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};
for i = 1:3
    nexttile(i)
    ax = gca;
    limits = [5 11];
    ticks = 5:1:11;
    ax.XLim = limits;
    ax.XTick = ticks;
    ax.XTickLabel = string(ticks);
    ax.FontSize = fs;
    ax.XLabel.FontSize = fs-1;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.025];
    ax.Box = 'off';
    text(0.02, 0.96, labels{i}, 'Units', 'normalized', 'FontSize', fs, 'FontWeight', 'bold')
    if i == 1
        ax.YLabel.String = 'Aitken N [cm^{-3}]';
    elseif i == 2
        ax.YLabel.String = 'Accumulation N [cm^{-3}]';
    elseif i == 3
        ax.YLabel.String = 'Supermicron N [cm^{-3}]';
    end
    ax.YLabel.FontSize = fs;
end

clim([0 25])
colormap(Colors)
c = colorbar;
c.Limits = [0 25];
c.Ticks = [0 5 10 15 20 25];
c.TickLabels = ["0" "5" "10" "15" "20" "25"];
c.Label.String = sprintf('Sea Surface Temperature\n (SST) [°C]');
c.Label.FontSize = fs;
c.TickDirection = 'out';
c.Layout.Tile = 'east';

t.XLabel.String = 'Subsurface Bubble Concentration (n) [cm^{-3}]';
t.XLabel.FontSize = fs;
t.YLabel.FontSize = fs;

set(gcf,'Position',[50 50 1150 300],'Color','w') % set standard figure size

%%% ------------------------------------------------------------------ %%%
% wieighted vs unweighted 

figure(5), clf
t = tiledlayout(2,4);

nexttile(1)
hold on
% combine data into table
yVar = "N_aitken";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x');
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2');
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.unweighted.(yVar).Lin_NRMSE, ...
    reg_metrics.unweighted.(yVar).Lin_p, ...
    reg_metrics.unweighted.(yVar).Lin_VarExplained, ...
    reg_metrics.unweighted.(yVar).Quad_NRMSE, ...
    reg_metrics.unweighted.(yVar).Quad_p, ...
    reg_metrics.unweighted.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.unweighted.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.unweighted.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end
ylabel('Non-weighted N [cm^{-3}]')

nexttile(2)
hold on
% combine data into table
yVar = "N_accumulation";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x');
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2');
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.unweighted.(yVar).Lin_NRMSE, ...
    reg_metrics.unweighted.(yVar).Lin_p, ...
    reg_metrics.unweighted.(yVar).Lin_VarExplained, ...
    reg_metrics.unweighted.(yVar).Quad_NRMSE, ...
    reg_metrics.unweighted.(yVar).Quad_p, ...
    reg_metrics.unweighted.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.unweighted.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.unweighted.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

nexttile(3)
hold on
% combine data into table
yVar = "N_supermicron";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x');
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2');
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.unweighted.(yVar).Lin_NRMSE, ...
    reg_metrics.unweighted.(yVar).Lin_p, ...
    reg_metrics.unweighted.(yVar).Lin_VarExplained, ...
    reg_metrics.unweighted.(yVar).Quad_NRMSE, ...
    reg_metrics.unweighted.(yVar).Quad_p, ...
    reg_metrics.unweighted.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.unweighted.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.unweighted.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

nexttile(4)
hold on
% combine data into table
yVar = "N";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x');
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2');
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.unweighted.(yVar).Lin_NRMSE, ...
    reg_metrics.unweighted.(yVar).Lin_p, ...
    reg_metrics.unweighted.(yVar).Lin_VarExplained, ...
    reg_metrics.unweighted.(yVar).Quad_NRMSE, ...
    reg_metrics.unweighted.(yVar).Quad_p, ...
    reg_metrics.unweighted.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.unweighted.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.unweighted.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

nexttile(5)
hold on
% combine data into table
yVar = "N_aitken";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.weighted.(yVar).Lin_NRMSE, ...
    reg_metrics.weighted.(yVar).Lin_p, ...
    reg_metrics.weighted.(yVar).Lin_VarExplained, ...
    reg_metrics.weighted.(yVar).Quad_NRMSE, ...
    reg_metrics.weighted.(yVar).Quad_p, ...
    reg_metrics.weighted.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.weighted.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.weighted.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end
ylabel('Weighted N [cm^{-3}]')

nexttile(6)
hold on
% combine data into table
yVar = "N_accumulation";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.weighted.(yVar).Lin_NRMSE, ...
    reg_metrics.weighted.(yVar).Lin_p, ...
    reg_metrics.weighted.(yVar).Lin_VarExplained, ...
    reg_metrics.weighted.(yVar).Quad_NRMSE, ...
    reg_metrics.weighted.(yVar).Quad_p, ...
    reg_metrics.weighted.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.weighted.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.weighted.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

nexttile(7)
hold on
% combine data into table
yVar = "N_supermicron";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.weighted.(yVar).Lin_NRMSE, ...
    reg_metrics.weighted.(yVar).Lin_p, ...
    reg_metrics.weighted.(yVar).Lin_VarExplained, ...
    reg_metrics.weighted.(yVar).Quad_NRMSE, ...
    reg_metrics.weighted.(yVar).Quad_p, ...
    reg_metrics.weighted.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.weighted.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.weighted.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

nexttile(8)
hold on
% combine data into table
yVar = "N";
for j = 1:size(data,2)
    tbl(j,:) = table(data{j}.T.mean, data{j}.T.mean.^2, ...
        data{j}.(yVar).mean, ...
        data{j}.(yVar).std, ...
        (1./data{j}.(yVar).std.^2));
end
tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

% linear model
mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
xFit = linspace(min(tbl.x), max(tbl.x), 100);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xFit;

% quadratic model
quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
xQuadFit = linspace(min(tbl.x), max(tbl.x), 100);
yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
    quad_mdl.Coefficients.Estimate(2) * xQuadFit + ...
    quad_mdl.Coefficients.Estimate(3) * xQuadFit.^2;

% run regression metrics, store result
[reg_metrics.weighted.(yVar).Lin_NRMSE, ...
    reg_metrics.weighted.(yVar).Lin_p, ...
    reg_metrics.weighted.(yVar).Lin_VarExplained, ...
    reg_metrics.weighted.(yVar).Quad_NRMSE, ...
    reg_metrics.weighted.(yVar).Quad_p, ...
    reg_metrics.weighted.(yVar).Quad_VarExplained] ...
    = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

% add linear and quadratic R^2
reg_metrics.weighted.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
reg_metrics.weighted.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

% Plot data and regression lines
hold on
plot(NaN, NaN, '--', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw)
plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw)

leg = legend;
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'off';

errorbar(tbl.x, tbl.y, tbl.std, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', lw-0.5)
plot(xFit, yFit, '--', 'Color', fColor, 'LineWidth', lw);
plot(xQuadFit, yQuadFit, ':', 'Color', fColor, 'LineWidth', lw);

for i = 1:size(data, 2)
    plot(tbl.x(i), tbl.y(i), LineStyle="none", Marker="o", MarkerFaceColor=Colors(i,:), MarkerEdgeColor="k")
end

labels = {'(a)\rm Aitken','(b)\rm Accumulation','(c)\rm Supermicron','(d)\rm Total','(e)','(f)', '(g)', '(h)'};

for i = 1:8
    nexttile(i)
    ax = gca;
    axList(i) = ax; 
    ax.XLabel.FontSize = fs;
    ax.YLabel.FontSize = fs;
    ax.FontSize = fs;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.025];
    text(0.02, 0.96, labels{i}, 'Units', 'normalized', 'FontSize', fs, 'FontWeight', 'bold')
end

YLimits = [
    0, 2500;
    0, 4000;
    0, 8;
    0, 6000];

for i = 1:4
    pair = [i, i+4];
    set(axList(pair), 'YLim', YLimits(i,:));
end

t.XLabel.String = 'Sea Surface Temperature (SST) [°C]';
t.XLabel.FontSize = fs;

set(gcf,'Position',[50 50 1150 650],'Color','w') % set standard figure size

%%% ------------------------------------------------------------------ %%%
%%% mean and averaged binned PNSDs

% load full data for final figure
cd('C:/Users/justi/OneDrive/Documents/GitHub/egusphere-2026-2142/')
full_smps_aps = load('soars_polar_smps_aps.mat').all_stats;
full_sm_ops = load('soars_polar_sm_ops.mat').all_stats;

% extract experiments 2 & 3 only
first_date = datetime(2024,02,11);
last_date = datetime(2024,03,02);

for i = 1:size(full_smps_aps, 2)
    idx(i) = full_smps_aps{i}.t.mean >= first_date & full_smps_aps{i}.t.mean <= last_date;
end

% remove experiments 1 & 2 from dataset prior to analysis
full_smps_aps = full_smps_aps(1, idx);
full_sm_ops = full_sm_ops(1, idx);

%%% Binned data
T_bins = 0:5:25;

for j = 1:length(T_bins)-1
    bin = [T_bins(j) T_bins(j+1)];
    for i = 1:size(full_smps_aps, 2)
        in_bin(i) = full_smps_aps{i}.T.mean > bin(1) & full_smps_aps{i}.T.mean <= bin(2);
    end
    new_data = full_smps_aps(1, in_bin);
    bin_data{j} = new_data;

    sm_ops_new = full_sm_ops(1, in_bin);
    sm_ops_bin{j} = sm_ops_new;
end
%%% ------------------------------------------------------------------ %%%
% SMPS-APS and SM-OPS Bin Mean and Binned Distributions
figure(6), clf

t = tiledlayout(2,5, "TileSpacing", "compact", "Padding", "compact");

for i = 1:size(data, 2)
    nexttile(i)
    hold on
    for j = 1:size(bin_data{i}, 2)
        plot(bin_data{i}{j}.D, bin_data{i}{j}.dN.mean, Color=Colors(i,:), LineWidth=lw)
    end
    plot(data{i}.D, data{i}.dN.mean, Color='k', LineWidth=lw)
end



for i = 1:size(smops, 2)
    nexttile(i+5)
    hold on
    for j = 1:size(sm_ops_bin{i}, 2)
        plot(sm_ops_bin{i}{j}.D, sm_ops_bin{i}{j}.dN.mean, Color=Colors(i,:), LineWidth=lw)
    end    
    plot(smops{i}.D, smops{i}.dN.mean, Color='k', LineWidth=lw)
end

labels = {'(a) \rmSMPS-APS','(b)','(c)','(d)', '(e)', ...
    '(f) \rmSM-OPS', '(g)', '(h)', '(i)', '(j)'};

for i = 1:10
    nexttile(i)
    text(0.02, 0.96, labels{i}, 'Units', 'normalized', 'FontSize', fs-2, 'FontWeight', 'bold')

    ax = gca;
    ax.XScale = 'log';  
    ax.FontSize = fs-2;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.025];
    ax.XTick = [10^-2 10^-1 10^0 10^1];
    ax.XLim = [10^-2 10^1];
    if i <= 5
        limits = [0 8000];
    else
        limits = [0 6000];
    end
    ax.YLim = limits;
    ax.YAxis.Exponent = 3;
    axis square
end

clim([0 25])
colormap(Colors)
c = colorbar;
c.Limits = [0 25];
c.Ticks = [0 5 10 15 20 25];
c.TickLabels = ["0" "5" "10" "15" "20" "25"];
c.Label.String = 'Sea Surface Temperature (SST) [°C]';
c.Label.FontSize = fs;
c.TickDirection = 'out';
c.Layout.Tile = 'east';

t.XLabel.String = 'Dry, Physical Diameter (D_{p, dry}) [\mum]';
t.XLabel.FontSize = fs;
t.YLabel.String = "PNSD (dN/dlogD_p) [cm^{-3}]";
t.YLabel.FontSize = fs;

set(gcf,'Position',[50 50 1150 550],'Color','w') % set standard figure size

%%% ------------------------------------------------------------------ %%%
% PDF of coldest and warmest SST N, S, and V distributions
figure(7), clf
t = tiledlayout(1,3);
nexttile(1)
hold on
legStr = [];
for i = [1 5]
    plot(data{i}.D, data{i}.dN.mean/trapz(log10(data{i}.D), data{i}.dN.mean), Color=Colors(i,:), LineWidth=lw)
    tempStr = sprintf(' %.1f', trapz(log10(data{i}.D), data{i}.dN.mean));
    legStr = [legStr; tempStr];
end
leg = legend;
leg.String = string(legStr);
leg.Box = 'off';
leg.Title.String = 'N [cm^{-3}]';
leg.Units = 'normalized';
leg.Position = [0.23 0.735 0.1 0.2];
leg.FontSize = fs-2;

nexttile(2)
hold on
legStr = [];
for i = [1 5]
    plot(data{i}.D, data{i}.dS.mean/trapz(log10(data{i}.D), data{i}.dS.mean), Color=Colors(i,:), LineWidth=lw)
    tempStr = sprintf(' %.1f', trapz(log10(data{i}.D), data{i}.dS.mean));
    legStr = [legStr; tempStr];
end
leg = legend;
leg.String = string(legStr);
leg.Box = 'off';
leg.Title.String = 'S [\mum^2 cm^{-3}]'; 
leg.Units = 'normalized';
leg.Position = [0.52 0.735 0.1 0.2];
leg.FontSize = fs-2;

nexttile(3)
hold on
legStr = [];
for i = [1 5]
    plot(data{i}.D, data{i}.dV.mean/trapz(log10(data{i}.D), data{i}.dV.mean), Color=Colors(i,:), LineWidth=lw)
    tempStr = sprintf(' %.1f', trapz(log10(data{i}.D), data{i}.dV.mean));
    legStr = [legStr; tempStr];
end
leg = legend;
leg.String = string(legStr);
leg.Box = 'off';
leg.Title.String = 'V [\mum^3 cm^{-3}]';  
leg.Units = 'normalized';
leg.Position = [0.81 0.735 0.1 0.2];
leg.FontSize = fs-2;

labels = {'(a)','(b)','(c)'};

for i = 1:3
    nexttile(i)
    ax = gca;
    ax.XScale = 'log';
    ax.XLim = [10^-2 10^1];
    ax.XTick = [10^-2 10^-1 10^0 10^1];
    ax.YLim = [0 2.5];
    ax.FontSize = fs;
    ax.YLabel.FontSize = fs-1;
    ax.XLabel.FontSize = fs-1;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.025];
    ax.Box = 'off';
    text(0.02, 0.96, labels{i}, 'Units', 'normalized', 'FontSize', fs, 'FontWeight', 'bold')
    axis square
end
t.XLabel.String = 'Dry, Physical Diameter (D_{p, dry}) [\mum]';
t.XLabel.FontSize = fs;
t.YLabel.String = 'PDF';
t.YLabel.FontSize = fs;

set(gcf,'Position',[50 50 1150 400],'Color','w') % set standard figure size

%% Save figures egusphere-2026-2142_Figure
cd('./figures/')
exportgraphics(figure(1), 'Figure3.png', ...
    'ContentType', 'image', 'Resolution', 300, 'Padding', 50)
exportgraphics(figure(2), 'FigureS7.png', ...
    'ContentType', 'image', 'Resolution', 300, 'Padding', 50)
exportgraphics(figure(3), 'Figure5.png', ...
    'ContentType', 'image', 'Resolution', 300, 'Padding', 50)
exportgraphics(figure(4), 'Figure6.png', ...
    'ContentType', 'image', 'Resolution', 300, 'Padding', 50)
exportgraphics(figure(5), 'FigureS4.png', ...
    'ContentType', 'image', 'Resolution', 300, 'Padding', 50)
exportgraphics(figure(6), 'FigureS6.png', ...
    'ContentType', 'image', 'Resolution', 300, 'Padding', 50)
exportgraphics(figure(7), 'FigureS8.png', ...
    'ContentType', 'image', 'Resolution', 300, 'Padding', 50)

%% PLC Figure
clear, clc, close all

% load in PLC for soars campaign 2-4
cd('C:/Users/justi/OneDrive/Documents/GitHub/egusphere-2026-2142/')
files = dir('*plc*.xlsx');
for i = 1:length(files)
    soars{i} = readtable(files(i).name);
    soars{i}.Properties.VariableNames = {'dp','efficiency','loss'};
end

% instrument index
spider_id = find(contains({files.name}, 'spider'));
smps_id = find(contains({files.name}, 'smps'));
aps_id = find(contains({files.name}, 'aps'));
ops_id = find(contains({files.name}, 'ops'));

%%% Figure Standardization
fs = 12; % set font size
lw = 1.5; % set line width
ms = 4; % set marker size

% instrument colors
sm_c = '#a6cee3';
smps_c = '#1f78b4';
ops_c = '#b2df8a';
aps_c = '#33a02c';

% Enter figure data
figure(1), clf

hold on
yline(100, LineStyle='--', LineWidth=lw-1, Color='k')

p1 = plot(soars{spider_id}.dp, soars{spider_id}.efficiency, ...
    LineStyle='-', LineWidth=lw, Color=sm_c);
p2 = plot(soars{smps_id}.dp, soars{smps_id}.efficiency, ...
    LineStyle='-', LineWidth=lw, Color=smps_c);
p3 = plot(soars{ops_id}.dp, soars{ops_id}.efficiency, ...
    LineStyle='-', LineWidth=lw, Color=ops_c);
p4 = plot(soars{aps_id}.dp, soars{aps_id}.efficiency, ...
    LineStyle='-', LineWidth=lw, Color=aps_c);

leg = legend([p1 p2 p3 p4]);
leg.String = {'SM', 'SMPS', 'OPS', 'APS'};
leg.AutoUpdate = 'off';
leg.Location = 'southwest';
leg.Box = 'on';

axis square


ax = gca;
ax.XScale = 'log';
ax.XLim = [10^-2 10^1];
ax.YLim = [0 120];
ax.FontSize = fs;
ax.TickDir = 'out';
ax.Box = 'off';
ax.TickLength = [0.02 0.025];
xline(ax.XLim(2))
yline(ax.YLim(2))
ax.XGrid = 'on';
ax.YGrid = 'on';


ax.XLabel.String = 'Dry, Physical Diameter (D_{p, dry}) [\mum]';
ax.XLabel.FontSize = fs;
ax.YLabel.String = 'Sampling Efficiency [%]';
ax.YLabel.FontSize = fs;
set(gcf,'Position',[50 50 500 500],'Color','w') % set standard figure size

%% Save figures 
cd('./figures/')
exportgraphics(figure(1), 'FigureS3.png', ...
    'ContentType', 'image', 'Resolution', 300, 'Padding', 50)

%% Hysteresis Plots
clear, clc, close all
cd('C:/Users/justi/OneDrive/Documents/GitHub/egusphere-2026-2142/')
data = load('soars_polar_final.mat').all_stats;
smops = load('soars_polar_smops_final.mat').all_stats;

% load full data 
full_smps_aps = load('soars_polar_smps_aps.mat').all_stats;
full_sm_ops = load('soars_polar_sm_ops.mat').all_stats;

% extract experiments 2 & 3 only
first_date = datetime(2024,02,11);
last_date = datetime(2024,03,02);

for i = 1:size(full_smps_aps, 2)
    idx(i) = full_smps_aps{i}.t.mean >= first_date & full_smps_aps{i}.t.mean <= last_date;
end

% remove experiments 1 & 2 from dataset prior to analysis
full_smps_aps = full_smps_aps(1, idx);
full_sm_ops = full_sm_ops(1, idx);

%%% Binned data
T_bins = 0:5:25;

for j = 1:length(T_bins)-1
    bin = [T_bins(j) T_bins(j+1)];
    for i = 1:size(full_smps_aps, 2)
        in_bin(i) = full_smps_aps{i}.T.mean > bin(1) & full_smps_aps{i}.T.mean <= bin(2);
    end
    new_data = full_smps_aps(1, in_bin);
    bin_data{j} = new_data;

    sm_ops_new = full_sm_ops(1, in_bin);
    sm_ops_bin{j} = sm_ops_new;
end

%%% Load Colormaps - From RJLIII
cd C:/Users/justi/OneDrive/Documents/0_UCSD/Research/MATLAB/Colormap/slanCM/
load('slanCM_Data.mat')
ColorMaps = slandarerCM;
clear('author','slandarerCM')
CM = flipud(ColorMaps(6).Colors{11});
for i = 1:5
    Colors(i,:) = CM(i*floor(256/5), 1:3);
end
Colors = flipud(Colors);

% reset file path
cd('C:/Users/justi/OneDrive/Documents/GitHub/egusphere-2026-2142/')
%%% Figure Standardization
fs = 14; % set font size
lw = 1.5; % set line width
ms = 8; % set marker size
fColor = "k";
labels = {"(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
    "(i)", "(j)", "(k)", "(l)"};

%% ------------------------------------------------------------------ %%%
dt_cutoff = datetime(2024,02,17);

yVars = ["N_aitken", "N_accumulation", "N_supermicron", "N"];
yStds = ["aitken_std", "accumulation_std", "supermicron_std", "N_std"];
cooling_rows = [];
warming_rows = [];

for i = 1:size(data, 2)
    for j = 1:size(bin_data{i}, 2)

        Tval   = bin_data{i}{j}.T.mean;
        tmean  = bin_data{i}{j}.t.mean;

        row = table( ...
            Tval, ...
            bin_data{i}{j}.N.mean, ...
            bin_data{i}{j}.N.std, ...
            bin_data{i}{j}.N_aitken.mean, ...
            bin_data{i}{j}.N_aitken.std, ...
            bin_data{i}{j}.N_accumulation.mean, ...
            bin_data{i}{j}.N_accumulation.std, ...
            bin_data{i}{j}.N_supermicron.mean, ...
            bin_data{i}{j}.N_supermicron.std, ...
            bin_data{i}{j}.S.mean, ...
            bin_data{i}{j}.S.std, ...
            bin_data{i}{j}.V.mean, ...
            bin_data{i}{j}.V.std, ...
            'VariableNames', ["T","N","N_std","N_aitken","aitken_std",...
            "N_accumulation","accumulation_std","N_supermicron",...
            "supermicron_std","S","S_std","V","V_std"] ...
            );

        if tmean <= dt_cutoff
            cooling_rows = [cooling_rows; row];
        else
            warming_rows = [warming_rows; row];
        end
    end
end

cooling = cooling_rows;
warming = warming_rows;

%%% plot
figure(1), clf
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.7], 'Color', 'w');

t = tiledlayout(4,3, "TileSpacing", "compact", "Padding", "compact");
t.YLabel.String = 'N [cm^{-3}]';
t.XLabel.String = 'Sea Surface Temperature (SST) [°C]';
t.XLabel.FontSize = fs;
t.YLabel.FontSize = fs;

coolTiles = [1 4 7 10];
warmTiles = [2 5 8 11];

for k = 1:numel(yVars)

    yVar = yVars(k);
    yStd = yStds(k);

    %%% -------------------------
    %  COOLING TILE
    %%% -------------------------
    nexttile(coolTiles(k))
    hold on

    % Plot cooling points
    for i = 1:size(data, 2)
        for j = 1:size(bin_data{i}, 2)

            Tval  = bin_data{i}{j}.T.mean;
            yval  = bin_data{i}{j}.(yVar).mean;
            tmean = bin_data{i}{j}.t.mean;

            if tmean <= dt_cutoff
                plot(Tval, yval, ...
                    Marker='v', MarkerFaceColor=Colors(i,:), ...
                    MarkerEdgeColor='k', LineStyle='none');
            end
        end
    end

    % Cooling regression
    tbl = table(cooling.T, cooling.T.^2, cooling.(yVar), ...
                cooling.(yStd), 1./cooling.(yStd).^2, ...
                'VariableNames', {'x','x2','y','std','w'});

    mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
    xFit = linspace(min(tbl.x), max(tbl.x), 200);
    yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*xFit;
    plot(xFit, yFit, 'k:', 'LineWidth', lw-0.5)

    quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
    yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
               quad_mdl.Coefficients.Estimate(2)*xFit + ...
               quad_mdl.Coefficients.Estimate(3)*xFit.^2;
    plot(xFit, yQuadFit, 'k-', 'LineWidth', lw-0.5)
    l1 = plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw-0.5);
    l2 = plot(NaN, NaN, '-', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw-0.5);

    leg = legend([l1 l2]);
    leg.AutoUpdate = 'off';
    leg.Location = 'best';
    leg.Box = 'off';
    
    %%% run regression metrics, store result
    [reg_metrics.cooling.(yVar).Lin_NRMSE, ...
        reg_metrics.cooling.(yVar).Lin_p, ...
        reg_metrics.cooling.(yVar).Lin_VarExplained, ...
        reg_metrics.cooling.(yVar).Quad_NRMSE, ...
        reg_metrics.cooling.(yVar).Quad_p, ...
        reg_metrics.cooling.(yVar).Quad_VarExplained] ...
        = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");
    
    % add linear and quadratic R^2
    reg_metrics.cooling.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
    reg_metrics.cooling.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

    %%% -------------------------
    %  WARMING TILE
    %%% -------------------------

    clear tbl
    nexttile(warmTiles(k))
    hold on

    % Plot warming points
    for i = 1:size(data, 2)
        for j = 1:size(bin_data{i}, 2)

            Tval  = bin_data{i}{j}.T.mean;
            yval  = bin_data{i}{j}.(yVar).mean;
            tmean = bin_data{i}{j}.t.mean;

            if tmean > dt_cutoff
                plot(Tval, yval, ...
                    Marker='^', MarkerFaceColor=Colors(i,:), ...
                    MarkerEdgeColor='k', LineStyle='none');
            end
        end
    end
    
    % Warming regression
    tbl = table(warming.T, warming.T.^2, warming.(yVar), ...
                warming.(yStd), 1./warming.(yStd).^2, ...
                'VariableNames', {'x','x2','y','std','w'});

    mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
    xFit = linspace(min(tbl.x), max(tbl.x), 200);
    yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*xFit;
    plot(xFit, yFit, 'k:', 'LineWidth', lw-0.5)

    quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
    yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
               quad_mdl.Coefficients.Estimate(2)*xFit + ...
               quad_mdl.Coefficients.Estimate(3)*xFit.^2;
    plot(xFit, yQuadFit, 'k-', 'LineWidth', lw-0.5)
    l1 = plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw-0.5);
    l2 = plot(NaN, NaN, '-', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw-0.5);

    leg = legend([l1 l2]);
    leg.AutoUpdate = 'off';
    leg.Location = 'best';
    leg.Box = 'off';

    %%% run regression metrics, store result
    [reg_metrics.warming.(yVar).Lin_NRMSE, ...
        reg_metrics.warming.(yVar).Lin_p, ...
        reg_metrics.warming.(yVar).Lin_VarExplained, ...
        reg_metrics.warming.(yVar).Quad_NRMSE, ...
        reg_metrics.warming.(yVar).Quad_p, ...
        reg_metrics.warming.(yVar).Quad_VarExplained] ...
        = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

    % add linear and quadratic R^2
    reg_metrics.warming.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
    reg_metrics.warming.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

    clear tbl
end

tOrder = [3 6 9 12];
for k = 1:numel(yVars)
    nexttile(tOrder(k))
    hold on
    yVar = yVars(k);
    yStd = yStds(k);
    % for each temperature bin
    for i = 1:size(data, 2)
        % for each sample within each bin
        for j = 1:size(bin_data{i}, 2)

            Tval  = bin_data{i}{j}.T.mean;
            yval  = bin_data{i}{j}.(yVar).mean;
            stdval = bin_data{i}{j}.(yVar).std;
            tmean = bin_data{i}{j}.t.mean;

            if tmean <= dt_cutoff
                mk = 'v';
            else
                mk = '^';
            end

            scatter(Tval, yval, ...
                36, Colors(i,:), ...
                'filled', ...
                'MarkerFaceAlpha', 0.35, ...
                'MarkerEdgeColor', fColor, ...
                'MarkerEdgeAlpha', 0.35, ...
                'Marker', mk);

            leg = legend([l1 l2]);
            leg.AutoUpdate = 'off';
            leg.Location = 'best';
            leg.Box = 'off';
        end

        errorbar(data{i}.T.mean, data{i}.(yVar).mean, ...
            data{i}.(yVar).std, data{i}.(yVar).std, ...
            Marker='o', MarkerFaceColor=Colors(i,:), ...
            MarkerEdgeColor='k', LineStyle='none', MarkerSize=ms, ...
            Color=fColor)

        % compile data into table for regressions
        tbl(i,:) = table(data{i}.T.mean, data{i}.T.mean.^2, ...
            data{i}.(yVar).mean, ...
            data{i}.(yVar).std, ...
            (1./data{i}.(yVar).std.^2));
        tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

    end

    %%% run regression metrics, store result
    [reg_metrics.combined.(yVar).Lin_NRMSE, ...
        reg_metrics.combined.(yVar).Lin_p, ...
        reg_metrics.combined.(yVar).Lin_VarExplained, ...
        reg_metrics.combined.(yVar).Quad_NRMSE, ...
        reg_metrics.combined.(yVar).Quad_p, ...
        reg_metrics.combined.(yVar).Quad_VarExplained] ...
        = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

    % add linear and quadratic R^2
    reg_metrics.combined.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
    reg_metrics.combined.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

    mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
    xFit = linspace(min(tbl.x), max(tbl.x), 200);
    yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*xFit;
    plot(xFit, yFit, 'k:', 'LineWidth', lw-0.5)

    quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
    yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
        quad_mdl.Coefficients.Estimate(2)*xFit + ...
        quad_mdl.Coefficients.Estimate(3)*xFit.^2;
    plot(xFit, yQuadFit, 'k-', 'LineWidth', lw-0.5)
    l1 = plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw-0.5);
    l2 = plot(NaN, NaN, '-', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw-0.5);

    leg = legend([l1 l2]);
    leg.AutoUpdate = 'off';
    leg.Location = 'best';
    leg.Box = 'off';
end

for i = 1:12
    nexttile(i)
    ax = gca;
    ax.FontSize = fs-2;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.025];
    ax.XLim = [0 25];
    text(0.02, 0.96, labels{i}, ...
        'Units', 'normalized', ...
        'FontSize', fs-2, ...
        'FontWeight', 'bold')

    if i == 1
        ax.YLabel.String = 'Aitken';
        ax.Title.String = 'Cooling';
    elseif i == 2
        ax.Title.String = 'Warming';
    elseif i == 3
        ax.Title.String = 'Combined';
    elseif i == 4
        ax.YLabel.String = 'Accumulation';
    elseif i == 7
        ax.YLabel.String = 'Supermicron';
    elseif i == 10
        ax.YLabel.String = 'Total';
    end
    if i < 7 || i > 9
        ax.YAxis.Exponent = 3;
    end
    if i <= 3
        ax.YLim = [0 2.5e3];
        ax.YTick = 0:1e3:2e3;
    elseif i == 4 || i == 5 || i == 6
        ax.YLim = [0 3.5e3];
        ax.YTick = 0:1e3:3e3;
    elseif i == 7 || i == 8 || i == 9
        ax.YLim = [0 8];
        ax.YTick = 0:2:8;
    elseif i > 9
        ax.YLim = [0 6e3];
        ax.YTick = 0:2e3:6e3;
    end
end

clim([0 25])
colormap(Colors)
c = colorbar;
c.Limits = [0 25];
c.Ticks = [0 5 10 15 20 25];
c.TickLabels = ["0" "5" "10" "15" "20" "25"];
c.Label.String = 'Sea Surface Temperature (SST) [°C]';
c.Label.FontSize = fs;
c.TickDirection = 'out';
c.Layout.Tile = 'east';

yVars = ["N","S","V"];
yStds = ["N_std","S_std","V_std"];

figure(2), clf
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.7], 'Color', 'w');

t = tiledlayout(3,3, "TileSpacing", "compact", "Padding", "compact");
t.XLabel.String = 'Sea Surface Temperature (SST) [°C]';
t.XLabel.FontSize = fs;

coolTiles = [1 4 7];
warmTiles = [2 5 8];

for k = 1:numel(yVars)

    yVar = yVars(k);
    yStd = yStds(k);

    %%% -------------------------
    %  COOLING TILE
    %%% -------------------------
    nexttile(coolTiles(k))
    hold on

    % Plot cooling points
    for i = 1:size(data, 2)
        for j = 1:size(bin_data{i}, 2)

            Tval  = bin_data{i}{j}.T.mean;
            yval  = bin_data{i}{j}.(yVar).mean;
            tmean = bin_data{i}{j}.t.mean;

            if tmean <= dt_cutoff
                plot(Tval, yval, ...
                    Marker='v', MarkerFaceColor=Colors(i,:), ...
                    MarkerEdgeColor='k', LineStyle='none');
            end
        end
    end

    % Cooling regression
    tbl = table(cooling.T, cooling.T.^2, cooling.(yVar), ...
                cooling.(yStd), 1./cooling.(yStd).^2, ...
                'VariableNames', {'x','x2','y','std','w'});

    mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
    xFit = linspace(min(tbl.x), max(tbl.x), 200);
    yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*xFit;
    plot(xFit, yFit, 'k:', 'LineWidth', lw-0.5)

    quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
    yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
               quad_mdl.Coefficients.Estimate(2)*xFit + ...
               quad_mdl.Coefficients.Estimate(3)*xFit.^2;
    plot(xFit, yQuadFit, 'k-', 'LineWidth', lw-0.5)
    l1 = plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw-0.5);
    l2 = plot(NaN, NaN, '-', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw-0.5);

    leg = legend([l1 l2]);
    leg.AutoUpdate = 'off';
    leg.Location = 'best';
    leg.Box = 'off';
    
    %%% run regression metrics, store result
    [reg_metrics.cooling.(yVar).Lin_NRMSE, ...
        reg_metrics.cooling.(yVar).Lin_p, ...
        reg_metrics.cooling.(yVar).Lin_VarExplained, ...
        reg_metrics.cooling.(yVar).Quad_NRMSE, ...
        reg_metrics.cooling.(yVar).Quad_p, ...
        reg_metrics.cooling.(yVar).Quad_VarExplained] ...
        = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");
    
    % add linear and quadratic R^2
    reg_metrics.cooling.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
    reg_metrics.cooling.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

    %%% -------------------------
    %  WARMING TILE
    %%% -------------------------

    clear tbl
    nexttile(warmTiles(k))
    hold on

    % Plot warming points
    for i = 1:size(data, 2)
        for j = 1:size(bin_data{i}, 2)

            Tval  = bin_data{i}{j}.T.mean;
            yval  = bin_data{i}{j}.(yVar).mean;
            tmean = bin_data{i}{j}.t.mean;

            if tmean > dt_cutoff
                plot(Tval, yval, ...
                    Marker='^', MarkerFaceColor=Colors(i,:), ...
                    MarkerEdgeColor='k', LineStyle='none');
            end
        end
    end
    
    % Warming regression
    tbl = table(warming.T, warming.T.^2, warming.(yVar), ...
                warming.(yStd), 1./warming.(yStd).^2, ...
                'VariableNames', {'x','x2','y','std','w'});

    mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
    xFit = linspace(min(tbl.x), max(tbl.x), 200);
    yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*xFit;
    plot(xFit, yFit, 'k:', 'LineWidth', lw-0.5)

    quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
    yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
               quad_mdl.Coefficients.Estimate(2)*xFit + ...
               quad_mdl.Coefficients.Estimate(3)*xFit.^2;
    plot(xFit, yQuadFit, 'k-', 'LineWidth', lw-0.5)
    l1 = plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw-0.5);
    l2 = plot(NaN, NaN, '-', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw-0.5);

    leg = legend([l1 l2]);
    leg.AutoUpdate = 'off';
    leg.Location = 'best';
    leg.Box = 'off';

    %%% run regression metrics, store result
    [reg_metrics.warming.(yVar).Lin_NRMSE, ...
        reg_metrics.warming.(yVar).Lin_p, ...
        reg_metrics.warming.(yVar).Lin_VarExplained, ...
        reg_metrics.warming.(yVar).Quad_NRMSE, ...
        reg_metrics.warming.(yVar).Quad_p, ...
        reg_metrics.warming.(yVar).Quad_VarExplained] ...
        = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

    % add linear and quadratic R^2
    reg_metrics.warming.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
    reg_metrics.warming.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

    clear tbl
end

tOrder = [3 6 9];
for k = 1:numel(yVars)
    nexttile(tOrder(k))
    hold on
    yVar = yVars(k);
    yStd = yStds(k);
    % for each temperature bin
    for i = 1:size(data, 2)
        % for each sample within each bin
        for j = 1:size(bin_data{i}, 2)

            Tval  = bin_data{i}{j}.T.mean;
            yval  = bin_data{i}{j}.(yVar).mean;
            stdval = bin_data{i}{j}.(yVar).std;
            tmean = bin_data{i}{j}.t.mean;

            if tmean <= dt_cutoff
                mk = 'v';
            else
                mk = '^';
            end

            scatter(Tval, yval, ...
                36, Colors(i,:), ...
                'filled', ...
                'MarkerFaceAlpha', 0.35, ...
                'MarkerEdgeColor', fColor, ...
                'MarkerEdgeAlpha', 0.35, ...
                'Marker', mk);

            leg = legend([l1 l2]);
            leg.AutoUpdate = 'off';
            leg.Location = 'best';
            leg.Box = 'off';
        end

        errorbar(data{i}.T.mean, data{i}.(yVar).mean, ...
            data{i}.(yVar).std, data{i}.(yVar).std, ...
            Marker='o', MarkerFaceColor=Colors(i,:), ...
            MarkerEdgeColor='k', LineStyle='none', MarkerSize=ms, ...
            Color=fColor)

        % compile data into table for regressions
        tbl(i,:) = table(data{i}.T.mean, data{i}.T.mean.^2, ...
            data{i}.(yVar).mean, ...
            data{i}.(yVar).std, ...
            (1./data{i}.(yVar).std.^2));
        tbl.Properties.VariableNames = {'x', 'x2', 'y', 'std', 'w'};

    end

    %%% run regression metrics, store result
    [reg_metrics.combined.(yVar).Lin_NRMSE, ...
        reg_metrics.combined.(yVar).Lin_p, ...
        reg_metrics.combined.(yVar).Lin_VarExplained, ...
        reg_metrics.combined.(yVar).Quad_NRMSE, ...
        reg_metrics.combined.(yVar).Quad_p, ...
        reg_metrics.combined.(yVar).Quad_VarExplained] ...
        = Regression_Metrics(tbl.x, tbl.y, mdl, quad_mdl, "fitlm");

    % add linear and quadratic R^2
    reg_metrics.combined.(yVar).Lin_R2 = mdl.Rsquared.Ordinary;
    reg_metrics.combined.(yVar).Quad_R2 = quad_mdl.Rsquared.Ordinary;

    mdl = fitlm(tbl, 'y ~ x', 'Weights', tbl.w);
    xFit = linspace(min(tbl.x), max(tbl.x), 200);
    yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*xFit;
    plot(xFit, yFit, 'k:', 'LineWidth', lw-0.5)

    quad_mdl = fitlm(tbl, 'y ~ x + x2', 'Weights', tbl.w);
    yQuadFit = quad_mdl.Coefficients.Estimate(1) + ...
        quad_mdl.Coefficients.Estimate(2)*xFit + ...
        quad_mdl.Coefficients.Estimate(3)*xFit.^2;
    plot(xFit, yQuadFit, 'k-', 'LineWidth', lw-0.5)
    l1 = plot(NaN, NaN, ':', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',mdl.Rsquared.Ordinary)], 'LineWidth', lw-0.5);
    l2 = plot(NaN, NaN, '-', 'Color', fColor, 'DisplayName', ['R^2 = ' sprintf('%.2f',quad_mdl.Rsquared.Ordinary)], 'LineWidth', lw-0.5);

    leg = legend([l1 l2]);
    leg.AutoUpdate = 'off';
    leg.Location = 'best';
    leg.Box = 'off';
end

for i = 1:9
    nexttile(i)
    ax = gca;
    ax.FontSize = fs-2;
    ax.TickDir = 'out';
    ax.XLim = [0 25];
    ax.TickLength = [0.02 0.025];
    text(0.02, 0.96, labels{i}, ...
        'Units', 'normalized', ...
        'FontSize', fs-2, ...
        'FontWeight', 'bold')

    if i == 1
        ax.YLabel.String = 'N [cm^{-3}]';
        ax.Title.String = 'Cooling';
    elseif i == 2
        ax.Title.String = 'Warming';
    elseif i == 3
        ax.Title.String = 'Combined';
    elseif i == 4
        ax.YLabel.String = 'S [\mum^2 cm^{-3}]';
    elseif i == 7
        ax.YLabel.String = 'V [\mum^3 cm^{-3}]';
    end

    if i <= 3
        ax.YLim = [0 6e3];
        ax.YTick = 0:2e3:6e3;
    elseif i == 4 || i == 5 || i == 6
        ax.YLim = [0 8e2];
        ax.YTick = 0:2e2:8e2;
    elseif i == 7 || i == 8 || i == 9
        ax.YLim = [0 50];
        ax.YTick = 0:10:50;
    end
end

clim([0 25])
colormap(Colors)
c = colorbar;
c.Limits = [0 25];
c.Ticks = [0 5 10 15 20 25];
c.TickLabels = ["0" "5" "10" "15" "20" "25"];
c.Label.String = 'Sea Surface Temperature (SST) [°C]';
c.Label.FontSize = fs;
c.TickDirection = 'out';
c.Layout.Tile = 'east';

%% Save figures egusphere-2026-2142_Figure
cd('./figures/')
exportgraphics(figure(1), 'FigureS10.png', ...
    'ContentType', 'image', 'Resolution', 300, 'Padding', 50)
exportgraphics(figure(2), 'FigureS11.png', ...
    'ContentType', 'image', 'Resolution', 300, 'Padding', 50)