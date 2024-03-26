function scheme = scheme_postC7(nu1)
% caluclates the POSTC7^1_2 sequence
%
% Output
%   seq.nu1   : rf amplitude (linear Hz)
%   seq.phi   : [list] rf phase of each pulse (degree)
%   seq.tau   : [list] duration of each pulse (sec.)
%
% Input
%   nu1   : rf amplitude (linear Hz)

% basic element
flip_element = [90,360,270]; % flips of basic element (degree)
phase_element = @(x) [x, mod(x+180,360),x]; 

% list of rf phases (degree)
ph = pi*[0,2/7,4/7,6/7,8/7,10/7,12/7]*180/pi; % (degree)
phases = cell2mat(arrayfun(phase_element,ph,'UniformOutput',false));

% pule duration array
flips = repmat(flip_element,1,7);
taus = flips/(360*nu1);                 

scheme.nu1 = nu1;
scheme.phi = phases;
scheme.tau = taus;

end