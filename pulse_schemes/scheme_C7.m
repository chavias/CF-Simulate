function scheme = scheme_C7(nu1)
% caluclates the C7 sequence
%
% Output
%   seq.nu1   : rf amplitude (linear Hz)
%   seq.phi   : [list] rf phase of each pulse (degree) DONE
%   seq.tau   : [list] duration of each pulse (sec.)
%
% Input
%   nu1   : rf amplitude (linear Hz)

% basic element
flip_element = [360,360]; % flips of basic element (degree)
phase_element = @(x) [x, mod(x+180,360)]; 

% list of rf phases (degree)
phases = pi*[0,2/7,4/7,6/7,8/7,10/7,12/7]*180/pi; % (degree)
phases = cell2mat(arrayfun(phase_element,phases,'UniformOutput',false));

% flip array
flips = repmat(flip_element,1,7);
taus = flips/(360*nu1);                 

scheme.nu1 = nu1;
scheme.phi = phases; 
scheme.tau = taus;


end