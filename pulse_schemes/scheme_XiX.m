function scheme = scheme_XiX(nu1)
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
flip_element = [20*180,20*180]; % flips of basic element (degree)
% list of rf phases (degree)
phi = [0,180]; % (degree)
% pule duration array
taus = flip_element./(360*nu1);                 

scheme.nu1 = nu1;
scheme.phi = phi;
scheme.tau = taus;

end