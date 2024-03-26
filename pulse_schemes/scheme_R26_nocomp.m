function scheme = scheme_R26_nocomp(nu1)
% caluclates the R26^11_4 sequence (N=26, n=4, nu=11)
%
% Output
%   seq.nu1   : rf amplitude (linear Hz)
%   seq.phi   : [list] rf phase of each pulse (degree)
%   seq.tau   : [list] duration of each pulse (sec.)
%
% Input
%   nu1   : rf amplitude (linear Hz)

nu = 11;
%phase = nu/26*pi
phase = nu/26*180; % (degree)

% basic element
flip_element = [180,180];  % (degree)
phase_element = @(x) mod([x, -x],360); % (degree)

% If the basic R element has a length of z pi rotations, the rf-
% field for the R sequence is w_1 = z N w_r / (2 n)

% sequence
flips = repmat(flip_element,1,26);  
phases = repmat(phase_element(phase),1,26); 
taus = flips/(360*nu1);

scheme.nu1 = nu1;
scheme.phi = phases;
scheme.tau = taus;

end