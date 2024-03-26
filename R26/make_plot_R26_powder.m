% create the graphics grid
clear; close all; clc 

%% load functions and data
addpath('./data/')
addpath('../../numerical_data/')
addpath('../../utilities/')


SetAllInterpreter2latex;
set(groot, 'DefaultLineLineWidth', 1);

load("R26_powder.mat")
time_resolution = 4.2735e-07;


f1 = figure('Name','R26 powder');
f1.Position(3:4) = [1400 800];
tiledlayout(2,3)

% Tile 1
load('R26O1r1.mat');
x_axis1 = O1.nur_list./O1.nu1;
x_axis = linspace(1e3,21e3,1001)./O1.nu1;

dataB = squeeze(data(:,2,index_T(O1.T,time_resolution)))./data(1,1,1);
nexttile
hold on
title(sprintf('$T= %.1f\\,$ms',round(O1.T*1e3,1)));
plot(x_axis1,O1.signalB2)
plot(x_axis,dataB)
hold off
ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(1,0)

% Tile 2
load('R26O1t1.mat');
dataB = squeeze(data(:,2,index_T(O1.T,time_resolution)))./data(1,1,1);
nexttile
hold on
title(sprintf('$T= %.1f\\,$ms',round(O1.T*1e3,1)));
plot(x_axis1,O1.signalB2)
plot(x_axis,dataB)
hold off
%ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(1,0)

% % Tile 3
load('R26O1t2.mat');
dataB = squeeze(data(:,2,index_T(O1.T,time_resolution)))./data(1,1,1);
nexttile
hold on
title(sprintf('$T= %.1f\\,$ms',round(O1.T*1e3,1)));
plot(x_axis1,O1.signalB2)
plot(x_axis,dataB)
hold off
%ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(1,1)

%% second-order

% Tile 4
load('R26O2r1.mat');
dataB = squeeze(data(:,2,index_T(O2.T,time_resolution)))./data(1,1,1);
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(O2.T*1e3,1)));
plot(x_axis1,O2.signalB2)
plot(x_axis,dataB)
hold off
ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(2,0)

% Tile 5
load('R26O2r2.mat');
dataB = squeeze(data(:,2,index_T(O2.T,time_resolution)))./data(1,1,1);
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(O2.T*1e3,1)));
plot(x_axis1,O2.signalB2)
plot(x_axis,dataB)
hold off
%ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(2,0)

% Tile 6
load('R26O2r3.mat');
dataB = squeeze(data(:,2,index_T(O2.T,time_resolution)))./data(1,1,1);
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(O2.T*1e3,1)));
plot(x_axis1,O2.signalB2)
plot(x_axis,dataB)
hold off
%ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(2,1)

%% enumerate

NumPlot(f1, {'(a)', '(b)', '(c)','(d)', '(e)', '(f)'}, 'VShift', 0, 'Direction', 'LeftRight', 'FontSize', 16)

%% export graphic
 
set(gcf, 'renderer', 'painters');
exportgraphics(gcf,'~/Documents/LaTeX/CF_effective/JCP/figures/R26_powder.pdf', ...
     'BackgroundColor','white','ContentType','vector');
 

%% utility functions

function [t1ms,t2ms,t3ms,x_axis] = load_data(data,time_resolution,order)
%function [t3ms,x_axis] = load_data(data,time_resolution,order)
% load analytical data
 
fprintf('Nr. 1 \n')
repetitions1 = 2;
t1ms = generate_C7_effective_MAS_powder(repetitions1,order);
t1ms.dataA = squeeze(data(:,1,index_T(t1ms.T,time_resolution)))./data(1,1,1);
t1ms.dataB = squeeze(data(:,2,index_T(t1ms.T,time_resolution)))./data(1,1,1);
 
fprintf('Nr. 2 \n')
repetitions2 = 6;
t2ms = generate_C7_effective_MAS_powder(repetitions2,order);
t2ms.dataA = squeeze(data(:,1,index_T(t2ms.T,time_resolution)))./data(1,1,1);
t2ms.dataB = squeeze(data(:,2,index_T(t2ms.T,time_resolution)))./data(1,1,1);
 
fprintf('Nr. 3 \n')
repetitions3 = 10;
t3ms = generate_C7_effective_MAS_powder(repetitions3,order);
t3ms.dataA = squeeze(data(:,1,index_T(t3ms.T,time_resolution)))./data(1,1,1);
t3ms.dataB = squeeze(data(:,2,index_T(t3ms.T,time_resolution)))./data(1,1,1);

x_axis = t3ms.nu1./t3ms.nur_list;
end

function index = index_T(T,resolution)
% calculates the index for a given time in data
index = round(T/resolution);
end

function setPlotParameterForOrder(order,leg)
% sets the parameter for plots
xlabel('$\nu_r/\nu_1$')
xlim([0.025,0.225])
%xlim([1/21,3/14])
%xticks([5.5,6,6.5,7,7.5,8,8.5])
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
