% First and second-order effective Hamiltonian for CW recoupling
% as a function of the rf amplitude nu1 for fixed MAS frequency.

function result = generate_CW_nu1_powder(T,order)

% set parameter for effective simulation
addpath('../../utilities/', "../../pulse_schemes/");

% experimental parameter
%T = 1.5e-3;                          % overall duration (sec.)

delta = -4.5e3;                       % dipolar coupling strength (linear Hz)
nucs = 0;                            % chemical-shift offset (linear Hz)
nur = 100e3;                         % MAS frequency
scheme.phi = 0;                      % scheme phase
scheme.tau = 100e-6;                 % leads to 1 rotation at 100kHz nu1 
scheme.nu1 = 10e3;
nu1_list = linspace(40e3,60e3,2001); % nu1 frequency list (linear Hz)  
% plot_scheme(scheme);               
thetam = acos(1/sqrt(3));            % magic angle (rad.)

% Fourier series coefficients parameter
npoints = 20;   % highest Fourier coefficients 
step = 1e-6;    % time resolution (sec.)
csflag = 1;     % include CS [0 False, 1 True]

% powder average parameter
count = 1:300; 
value1= 300;
value2= 37;
value3= 61;
beta  = pi * count/value1;
alpha = 2*pi * rem((value2*count),value1)/value1;
gamma = 2*pi * rem((value3*count),value1)/value1;

%% definitions and allocations

%tau_m = sum(scheme.tau); % modulation time (sec.)

% definition of the spin operators 
[Ix,Iy,Iz, ...
 Ixx,Ixy,Ixz, ...
 Iyx,Iyy,Iyz, ...
 Izx,Izy,Izz, ...
 Ix1,Iy1,Iz1, ...
 I1x,I1y,I1z, ...
 Iv] = spin_operators;

% calculate the Fourier coefficients
% pulse scheme and cs offset is constant
% nu_r changes
[coeff,~,~,~] = sequence_get_coeff(...
    scheme.tau, ...   % duration of each pulse
    scheme.phi, ...   % rf phase of each pulse
    scheme.nu1, ...   % rf amplitude of pulse scheme (Hz linear)
    nucs, ...         % cs offset (Hz linear)
    step, ...         % time resolution (sec.)
    npoints, ...      % number of coefficients
    csflag);
% plot_coefficients(coeff,2)
% datax2 is normalized: 1/sqrt(6) (2zz-xx-yy)
[~, datax2, ~, ~, ~] = general_calc_twospin(coeff);

% data allocation
signal2A1=zeros(size(nu1_list));
signal2A2=zeros(size(nu1_list));
signal2B1=zeros(size(nu1_list));
signal2B2=zeros(size(nu1_list));

% waitbar
msg = sprintf("Up to second-order with %d powder points", count(end));
f2 = waitbar(0,msg,'Name','Second-order',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f2,'canceling',0);

for nu1_index=1:length(nu1_list) % nur loop
    % waitbar
    if getappdata(f2,'canceling')
        break
    end
    waitbar(nu1_index/numel(nu1_list),f2)

    for powder_index=1:length(beta) % powder loop
        scale=sin(beta(powder_index));
        % spatial part: coupling coefficients
        wn = zeros(1,5);
        for n=-2:2
            wn(n+3) = dlmn(2,n,0,thetam)*dlmn(2,0,n,beta(powder_index)) ...
                *exp(-1i*n*gamma(powder_index))*sqrt(3/2)*delta;
        end
        % possible optimization: could be outside the loop
        % spin part: two-spin coefficients * spin operator
        Hn = zeros(4,4,5);
        for m=-2:2:2
            for s=1:9 % spin operators
                Hn(:,:,m+3) = Hn(:,:,m+3) + datax2(npoints+m,3,s)*Iv(:,:,s);
            end
        end
        % full Hamiltonian
        HamA = zeros(4,4);
        for m=-2:2:2 % indexing wm
            for n=-2:1:2
                % effective frequency is not present
                nuA = n*nur + m*nu1_list(nu1_index); % linear in Hz
                HamA = HamA + Hn(:,:,m+3)*wn(n+3)*calc_h1(nuA,T);
            end
        end
        if order == 2
            % second-order effective Hamiltonian
            for m1 = -2:2:2 % coefficients A index
                for m2 = -2:2:2 % coefficients B index
                    comm = commutator(Hn(:,:,m1+3),Hn(:,:,m2+3));
                    for n2=-2:1:2 % dipolar coupling B
                        for n1=-2:1:2 % dipolar coupling A
                            % resonance condition A (n1,n2)
                            nuA = m1*nu1_list(nu1_index)+n1*nur; % Hz linear
                            % resonance condition B (m1,m2)
                            nuB = m2*nu1_list(nu1_index)+n2*nur; % Hz linear
                            % Commutator of the Fourier coefficients and coupling strength
                            Ham2 = -1/2*wn(n1+3)*wn(n2+3)*comm;
                            % second-order basis function
                            h2 = calc_h2(nuA,nuB,T);
                            Ham2 = Ham2 .* h2;
                            HamA = HamA + Ham2;
                        end
                    end
                end
            end
        end
    % the first and second-order Hamiltonian
    UA  = expm(-1i*2*pi*(HamA)*T);
    UAi = expm(+1i*2*pi*(HamA)*T);

    % detection
    sigmaA = I1x;
    det1   = I1x;
    det2   = Ix1;
    sigmaB = UA*sigmaA*UAi;

    signal2A1(nu1_index) = signal2A1(nu1_index)+trace(sigmaA*det1)*scale;
    signal2A2(nu1_index) = signal2A2(nu1_index)+trace(sigmaA*det2)*scale;
    signal2B1(nu1_index) = signal2B1(nu1_index)+trace(sigmaB*det1)*scale;
    signal2B2(nu1_index) = signal2B2(nu1_index)+trace(sigmaB*det2)*scale;
    
    end % powder loop
end % w1 loop
delete(f2); % remove waitbar

% normalize intensities with initial spin 1 intenstiy
result.signal2A2 = real(signal2A2/signal2A1(1));
result.signal2B2 = real(signal2B2/signal2A1(1));
result.signal2B1 = real(signal2B1/signal2A1(1));
result.signal2A1 = real(signal2A1/signal2A1(1));

result.nu1_list = nu1_list;
result.nur = nur;

end

%% functions

function h = calc_h1(nuA,T)
% Basis function of the first-order effective Hamiltonian
% nuA: resonance condition A (Hz linear)

% sinc(x) = sin(pi*x)/(pi*x)
h=sinc(2*nuA*T/2);

end

function h = calc_h2(nuA,nuB,T)
% Basis function of the second-order effective Hamiltonian
% Expressions are derived in second-order-hamiltonian.nb
% nuA:  resonance condition A (Hz linear)
% nuB:  resonance condition B (Hz linear)
% T:   duration (sec.)
wA = nuA*2*pi;
wB = nuB*2*pi;
threshold = 1e-7;
if abs(wA) < threshold
    h = wB.^(-2)*(wB*T*cos((1/2)*wB*T)+(-2)*sin((1/2)*wB*T));
elseif abs(wB) < threshold
    h = (-1)*wA.^(-2)*(wA*T*cos((1/2)*wA*T)+(-2)*sin((1/2)*wA*T));
elseif abs(wA + wB) < threshold
    h = wA.^(-2)*(wA*T+(-1)*sin(wA*T));
else
    h = 2/(wA*wB*(wA+wB))*(wB*cos(1/2*wB*T)* ...
        sin(1/2*wA*T)-wA*cos(1/2*wA*T)*sin(1/2*wB*T));
end
h = 2*pi*h/T;
if isnan(h)
    h = 0;
end
end
