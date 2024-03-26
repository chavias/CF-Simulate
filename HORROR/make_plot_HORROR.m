% create the graphics grid
clear; close all

addpath('./data/')
addpath('../../utilities/')

%% load numerical data
time_resolution = 2e-7;

data_file = 'HORROR_powder.mat';
numerical = load(data_file,'data');

%% load analytical data

% first-order data 
%[O1_t1,O1_t2,O1_t3,~]  = load_data(numerical.data,time_resolution,1);
% second-order data 
%[O2_t1,O2_t2,O2_t3,x_axis]  = load_data(numerical.data,time_resolution,2);
%save('analytical_HORROR_powder.mat','O1_t1','O1_t2','O1_t3','O2_t1','O2_t2','O2_t3','x_axis');


load('analytical_HORROR_powder.mat');

%% make figure

SetAllInterpreter2latex;
set(groot, 'DefaultLineLineWidth', 1);

T1 = 0.25e-3;
T2 = 0.5e-3;
T3 = 1e-3;

f1 = figure('Name','HORROR');
f1.Position(3:4) = [1400 800];
tiledlayout(2,3)

% Tile 1
nexttile
hold on
title(sprintf('$T= %.1f\\,$ms',round(T1*1e3,1)));
p1 = plot(x_axis,O1_t1.signal2B2);
p2 = plot(x_axis,O1_t1.dataB);
hold off
ylabel('$\langle \mathrm{I}_{2x} \rangle$')
setPlotParameterForOrder(1,0)
% Tile 2
nexttile
hold on
title(sprintf('$T= %.1f\\,$ms',round(T2*1e3,1)));
p = plot(x_axis,O1_t2.signal2B2);
p = plot(x_axis,O1_t2.dataB);
hold off
setPlotParameterForOrder(1,0)
% Tile 3
nexttile
title(sprintf('$T= %.1f\\,$ms',round(T3*1e3,1)));
hold on
p = plot(x_axis,O1_t3.signal2B2);
p = plot(x_axis,O1_t3.dataB);
hold off
setPlotParameterForOrder(1,1)

% Tile 4
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(T1*1e3,1)));
p = plot(x_axis,O2_t1.signal2B2);
p = plot(x_axis,O2_t1.dataB);
hold off
ylabel('$\langle \mathrm{I}_{2x} \rangle$')
setPlotParameterForOrder(2,0)
% Tile 5
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(T2*1e3,1)));
p = plot(x_axis,O2_t2.signal2B2);
p = plot(x_axis,O2_t2.dataB);
hold off
setPlotParameterForOrder(2,0)
% Tile 6
nexttile
%title(sprintf('$T= %.1f\\,$ms',round(T3*1e3,1)));
hold on
p = plot(x_axis,O2_t3.signal2B2);
p = plot(x_axis,O2_t3.dataB);
hold off
setPlotParameterForOrder(2,1)

% save figure
% orient(f1,'landscape')
% figure_name = sprintf('../figures/figure_HORROR');
% savefig(strcat(figure_name,'.fig')) 
% print(strcat(figure_name,'.pdf'), '-bestfit','-dpdf')  

%% numerate

NumPlot(f1, {'(a)', '(b)', '(c)','(d)', '(e)', '(f)'}, 'VShift', 0, 'Direction', 'LeftRight', 'FontSize', 16)

%% export

set(gcf, 'renderer', 'painters');
exportgraphics(gcf,'~/Documents/LaTeX/CF_effective/figures/HORROR_powder.pdf', ...
    'BackgroundColor','white','ContentType','vector');

%% utility functions

function [t1ms,t2ms,t3ms,x_axis] = load_data(data,time_resolution,order)

T1 = 0.25e-3;
t1ms = generate_CW_nu1_powder(T1,order);
t1ms.dataA = squeeze(data(:,1,index_T(T1,time_resolution)))./data(1,1,1);
t1ms.dataB = squeeze(data(:,2,index_T(T1,time_resolution)))./data(1,1,1);

T2 = 0.5e-3;
t2ms = generate_CW_nu1_powder(T2,order);
t2ms.dataA = squeeze(data(:,1,index_T(T2,time_resolution)))./data(1,1,1);
t2ms.dataB = squeeze(data(:,2,index_T(T2,time_resolution)))./data(1,1,1);

T3 = 1e-3;
t3ms = generate_CW_nu1_powder(T3,order);
t3ms.dataA = squeeze(data(:,1,index_T(T3,time_resolution)))./data(1,1,1);
t3ms.dataB = squeeze(data(:,2,index_T(T3,time_resolution)))./data(1,1,1);

x_axis = t1ms.nu1_list/t1ms.nur;

end

function index = index_T(T,resolution)
% calculates the index for a given time in data
index = round(T/resolution);
end

function setPlotParameterForOrder(order,leg)
% sets the parameter for plots
xlabel('$\nu_1/\nu_r$')
xlim([0.4,0.6])
%xticks([0.47,0.48 0.49 0.5 0.51 0.52,0.53])
if leg == 1
    switch true
        case order == 1
            legend('$\bar{\mathrm{H}}_{\mathrm{eff}}^{(1)}$','exact', 'Location','SouthEast')
        case order == 2
            legend('$\bar{\mathrm{H}}_{\mathrm{eff}}^{(1)}+\bar{\mathrm{H}}_{\mathrm{eff}}^{(2)}$','exact' ...
                ,'Location','SouthEast')
    end
end
set(gca,'FontSize',16)
ylim([-1,0])
grid on
set(gca,'LineWidth',1);
box on
end

