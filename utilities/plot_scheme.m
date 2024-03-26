function plot_scheme(seq, nu_r)
% plots the amplitude phase of the sequence

nu1= real(seq.nu1); % rf amplitude (linear)
phi = real(seq.phi);
tau = seq.tau;

% create amplitude, phase, and time vector
time_resolution = 1e-6;
tEnd = 0;
total_steps = sum(round(tau ./ time_resolution));
nu1_vector = zeros(1, total_steps);
phi_vector = zeros(1, total_steps);
time_vector = zeros(1, total_steps);
currentIndex = 1;
for i = 1:numel(phi)
    n_steps = round(tau(i) / time_resolution);
    nu1_pulse = repmat(nu1, 1, n_steps);
    phi_pulse = repmat(phi(i), 1, n_steps);
    t1 = tEnd;
    tEnd = t1 + tau(i);
    time_pulse = linspace(t1, tEnd, n_steps);
    endIndex = currentIndex + n_steps - 1;
    nu1_vector(currentIndex:endIndex) = nu1_pulse;
    phi_vector(currentIndex:endIndex) = phi_pulse;
    time_vector(currentIndex:endIndex) = time_pulse;
    currentIndex = endIndex + 1;
end

% convert to rf amplitude to kHz
nu1_vector_kHz = nu1_vector./1e3;
time_vector_ms = time_vector.*1e3;

% create figure
figure;
SetAllInterpreter2latex
% plot amplitude
subplot(2,1,1)
plot(time_vector_ms,nu1_vector_kHz);
% plot MAS frequency if given
if exist('nu_r','var')
    % second parameter does not exist, so default it to something
    xline(1/(nu_r/1e3), 'r-');
end
set(gca,'XLim',[0,time_vector_ms(end)]);
hXLabel = xlabel('$t$/ms');
hYLabel = ylabel('$\nu_1$ /kHz');
SetPlotStyle(hXLabel,hYLabel)

% plot phase
subplot(2,1,2);
plot(time_vector_ms,phi_vector);
if exist('nu_r','var')
    % second parameter does not exist, so default it to something
    xline(1/(nu_r/1e3), 'r-','\tau_r');
end
ylim_low  = round(min(phi_vector./10))*10-30;
ylim_high = round(max(phi_vector./10))*10+30;
set(gca,'XLim',[0,time_vector_ms(end)]);
set(gca,'YLim',[ylim_low, ylim_high]);
%yticks(linspace(-180,180,37));
hXLabel = xlabel('$t$/ms');
hYLabel = ylabel('$\phi/^\circ$');
%sgtitle(name)
SetPlotStyle(hXLabel,hYLabel)
end

function SetAllInterpreter2latex
% sets all interpreter to latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
end

function SetPlotStyle(hXLabel, hYLabel)
set([hXLabel, hYLabel], 'FontName', 'AvantGarde','FontSize', 14)
set(gcf,'color','w');
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
    'XColor', [0 0 0], 'YColor', [0 0 0],'fontsize',12,'FontWeight','normal', ...
    'LineWidth', 0.2)
end