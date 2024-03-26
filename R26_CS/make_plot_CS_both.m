% create the graphics grid
clear; close all; clc

%% load path

addpath('../../utilities/')
addpath('../../pulse_schemes/')

%% plot

f1 = figure('Name','R26 chemical shift');
f1.Position(3:4) = [1400 800];
tiledlayout(2,3)


%% figure

SetAllInterpreter2latex;
set(groot, 'DefaultLineLineWidth', 1);

time_resolution = 4.2735e-07;

rep = [2,6,10];

%% uncompensated

addpath('../R26_uncomp/')
load('R26_uncomp_powder_CS.mat')

load('N26O2t1CS.mat');
x_axis = O2.nucs_list*1e-3./70;
dataB = squeeze(data(:,2,index_T(O2.T,time_resolution)))./data(1,1,1);
nexttile
hold on
title(sprintf('$T= %.1f\\,$ms',round(O2.T*1e3,1)));
plot(x_axis,O2.signalB2)
plot(x_axis,dataB)
hold off
ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(2,0)

% Tile 5
load('N26O2t2CS.mat');
dataB = squeeze(data(:,2,index_T(O2.T,time_resolution)))./data(1,1,1);
nexttile
hold on
title(sprintf('$T= %.1f\\,$ms',round(O2.T*1e3,1)));
plot(x_axis,O2.signalB2)
plot(x_axis,dataB)
hold off
%ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(2,0)

% Tile 6
load('N26O2t3CS.mat');
dataB = squeeze(data(:,2,index_T(O2.T,time_resolution)))./data(1,1,1);
nexttile
hold on
title(sprintf('$T= %.1f\\,$ms',round(O2.T*1e3,1)));
plot(x_axis,O2.signalB2)
plot(x_axis,dataB)
hold off
ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(2,1)


%% compensated

% numerical
load("R26_powder_CS.mat")

% Tile 1
load('R26O1r1CS.mat');
x_axis = O2.nucs_list*1e-3./70;
%x_axis = O1.nucs_list*1e-3;
dataB = squeeze(data(:,2,index_T(O1.T,time_resolution)))./data(1,1,1);
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(O1.T*1e3,1)));
plot(x_axis,O1.signalB2)
plot(x_axis,dataB)
hold off
ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(2,0)

% Tile 2
load('R26O1t1CS.mat');
dataB = squeeze(data(:,2,index_T(O1.T,time_resolution)))./data(1,1,1);
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(O1.T*1e3,1)));
plot(x_axis,O1.signalB2)
plot(x_axis,dataB)
hold off
%ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(2,0)

% % Tile 3
load('R26O1t2CS.mat');
dataB = squeeze(data(:,2,index_T(O1.T,time_resolution)))./data(1,1,1);
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(O1.T*1e3,1)));
plot(x_axis,O1.signalB2)
plot(x_axis,dataB)
hold off
%ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(2,1)

%% enumerate

NumPlot(f1, {'(a)', '(b)', '(c)','(d)', '(e)', '(f)'},'Hshift', -0.003 ,'VShift', 0, 'Direction', 'LeftRight', 'FontSize', 16)


%% export figure

set(gcf, 'renderer', 'painters');
exportgraphics(gcf,'~/Documents/LaTeX/CF_effective/JCP/figures/R26_powder_CS_both.pdf', ...
    'BackgroundColor','white','ContentType','vector');


%% functions

function index = index_T(T,resolution)
% calculates the index for a given time in data
index = round(T/resolution);
end

function setPlotParameterForOrder(order,leg)
% sets the parameter for plots
%xlabel('$\nu_\mathrm{cs}/\omega_{\mathrm{I}}^{(0)}$kHz')
xlabel('$\nu_{\mathrm{I}}^{(0)}/\nu_1$')
%xlim([-1,1])
if leg == 1
    switch true
        case order == 1
            legend('$\bar{\mathrm{H}}_{\mathrm{eff}}^{(1)}$', ...
                'exact', 'Location','SouthEast')
        case order == 2
            legend(['$\bar{\mathrm{H}}_{\mathrm{eff}}^{(1)}' ...
                '+\bar{\mathrm{H}}_{\mathrm{eff}}^{(2)}$'],'exact' ...
                ,'Location','SouthEast')
    end
end
set(gca,'FontSize',16)
ylim([-1,0])
grid on
box on
end

function SetAllInterpreter2latex
% Sets all interpreter to latex
%
% Author: Matias Chavez
% Date: 12.04.22
    list_factory = fieldnames(get(groot,'factory'));
    index_interpreter = find(contains(list_factory,'Interpreter'));
    for i = 1:length(index_interpreter)
        default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
        set(groot,default_name,'latex');
    end
end


