function plot_coefficients(coeff,xyz)
%PLOT_COEFFICIENTS plots the single-spin coefficients
% Input: 
%       coeff : fourier coeff. [2*n-1, 3, 9]
%       xyz (optional) : effective field resonance

% determine which figure(s) to create
if exist('xyz','var')
   i_index = determine_which_plot(xyz);
else
   % plot all
   i_index = 1:3;
end
% x axis
x_coeff = calc_coeff_number(coeff);
% make figures
for i=i_index
    plot_names = {'-1','0','1'};
    figure('Name',sprintf('Plot %s',plot_names{i}));
    plot_subplot(x_coeff,coeff,i)
    sgtitle(sprintf('%d $\\nu_{\\mathrm{eff}}$',i-2))
end
end

%% utilities

function x_coeff = calc_coeff_number(coeff)
    % determine number of coefficients
    scale=numel(coeff(:,1,1))-1;
    x_coeff = (-scale/2):(scale/2);
end

function plot_subplot(x_coeff,coeff,i)
    for j=1:9
        subplot(3,3,j)
%        scatter(x_coeff,squeeze(real(coeff(:,i,j))), ...
%            'filled','o', ...
%            'SizeData', 10)
        plot(x_coeff,squeeze(real(coeff(:,i,j))))
        xlim([x_coeff(1),x_coeff(end)])
        ylim([-0.2,0.4])
    end
end

function i_index=determine_which_plot(xyz)
switch xyz
    case {'x' , -1}
        i_index = 1;
    case {'y' , 0}
        i_index = 2;
    case {'z' , 1}
        i_index = 3;
    otherwise 
        i_index = 1:3;
end
end
