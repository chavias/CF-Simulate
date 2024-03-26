% create the graphics grid
clear; close all; clc 

%% load functions

addpath('./data/')
addpath('../../utilities/')

reps = [2,6,10];

load("C7_powder.mat")
time_resolution = 1.98413e-07;

rep = [2,6,10];

SetAllInterpreter2latex;
set(groot, 'DefaultLineLineWidth', 1);

f1 = figure('Name','postC7 chemical shift');
f1.Position(3:4) = [1400 800];
tiledlayout(2,3)

% Tile 1
load('C7O1t1.mat');
x_axis =O1.nur_list./O1.nu1;
dataB = squeeze(data(:,2,index_T(O1.T,time_resolution)))./data(1,1,1);
nexttile
hold on
title(sprintf('$T= %.1f\\,$ms',round(O1.T*1e3,1)));
plot(x_axis,O1.signalB2)
plot(x_axis,dataB)
hold off
%ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(1,0)

% Tile 2
load('C7O1t2.mat');
dataB = squeeze(data(:,2,index_T(O1.T,time_resolution)))./data(1,1,1);
nexttile
hold on
title(sprintf('$T= %.1f\\,$ms',round(O1.T*1e3,1)));
plot(x_axis,O1.signalB2)
plot(x_axis,dataB)
hold off
%ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(1,0)

% % Tile 3
load('C7O1t3.mat');
dataB = squeeze(data(:,2,index_T(O1.T,time_resolution)))./data(1,1,1);
nexttile
hold on
title(sprintf('$T= %.1f\\,$ms',round(O1.T*1e3,1)));
plot(x_axis,O1.signalB2)
plot(x_axis,dataB)
hold off
%ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(1,1)

% Tile 4
load('C7O2t1.mat');
dataB = squeeze(data(:,2,index_T(O2.T,time_resolution)))./data(1,1,1);
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(O2.T*1e3,1)));
plot(x_axis,O2.signalB2)
plot(x_axis,dataB)
hold off
ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(2,0)

% Tile 5
load('C7O2t2.mat');
dataB = squeeze(data(:,2,index_T(O2.T,time_resolution)))./data(1,1,1);
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(O2.T*1e3,1)));
plot(x_axis,O2.signalB2)
plot(x_axis,dataB)
hold off
%ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(2,0)

% Tile 6
load('C7O2t3.mat');
dataB = squeeze(data(:,2,index_T(O2.T,time_resolution)))./data(1,1,1);
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(O2.T*1e3,1)));
plot(x_axis,O2.signalB2)
plot(x_axis,dataB)
hold off
%ylabel('$\langle \mathrm{I}_{2z} \rangle$')
setPlotParameterForOrder(2,1)

%% numerate

%addpath('/home/mach/Documents/MATLAB/C7_effective/utilities')
NumPlot(f1, {'(a)', '(b)', '(c)','(d)', '(e)', '(f)'}, 'VShift', 0, 'Direction', 'LeftRight', 'FontSize', 16)

%% export
set(gcf, 'renderer', 'painters');
exportgraphics(gcf,'~/Documents/LaTeX/CF_effective/figures/C7_powder.pdf', ...
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
