f1 = figure('Name','rotary resonance');
f1.Position(3:4) = [1400 800];
tiledlayout(2,3)

% Tile 1
% nexttile
% hold on
% title(sprintf('$T= %.1f\\,$ms',round(T.T1*1e3,1)));
% plot(x_axis,O1_t1.signal2B1)
% plot(x_axis,O1_t1.dataA)
% hold off
% ylabel('$\langle \mathrm{I}_{2x} \rangle$')
% setPlotParameterForOrder(1,0)
% % Tile 2
% nexttile
% hold on
% title(sprintf('$T= %.1f\\,$ms',round(T.T2*1e3,1)));
% plot(x_axis,O1_t2.signal2B1)
% plot(x_axis,O1_t2.dataA)
% hold off
% setPlotParameterForOrder(1,0)
% % Tile 3
% nexttile
% title(sprintf('$T= %.1f\\,$ms',round(T.T3*1e3,1)));
% hold on
% plot(x_axis,O1_t3.signal2B1)
% plot(x_axis,O1_t3.dataA)
% hold off
% setPlotParameterForOrder(1,1)

% Tile 4
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(T1*1e3,1)));
plot(x_axis,O2_t1.signal2B1)
plot(x_axis,O2_t1.dataA)
hold off
ylabel('$\langle \mathrm{I}_{2x} \rangle$')
setPlotParameterForOrder(2,0)
% Tile 5
nexttile
hold on
%title(sprintf('$T= %.1f\\,$ms',round(T2*1e3,1)));
plot(x_axis,O2_t2.signal2B1)
plot(x_axis,O2_t2.dataA)
hold off
setPlotParameterForOrder(2,0)
% Tile 6
nexttile
%title(sprintf('$T= %.1f\\,$ms',round(T3*1e3,1)));
hold on
plot(x_axis,O2_t3.signal2B1)
plot(x_axis,O2_t3.dataA)
hold off
setPlotParameterForOrder(2,1)

% orient(f1,'landscape')
% figure_name = sprintf('../figures/figure_rotary_resonance');
% savefig(strcat(figure_name,'.fig')) 
% print(strcat(figure_name,'.pdf'), '-bestfit','-dpdf') 

%% function

function index = index_T(T,resolution)
% calculates the index for a given time in data
index = round(T/resolution);
end

function setPlotParameterForOrder(order,leg)
% sets the parameter for plots
xlabel('$\nu_1/\nu_r$')
%xlim([0.98,1.03])
%xticks([0.95,0.96,0.97,0.98,0.99,1,1.01,1.02,1.03,1.04,1.05])
%xticks([0.97,0.98,0.99,1,1.01,1.02,1.03])
if leg == 1
    switch true
        case order == 1
            legend('$\bar{\mathrm{H}}_{\mathrm{eff}}^{(1)}$','exact', 'Location','SouthEast')
        case order == 2
            legend('$\bar{\mathrm{H}}_{\mathrm{eff}}^{(1)}+\bar{\mathrm{H}}_{\mathrm{eff}}^{(2)}$','exact' ...
                ,'Location','SouthEast')
    end
end
set(gca,'FontSize',14)
ylim([-1,1])
grid on
box on
end