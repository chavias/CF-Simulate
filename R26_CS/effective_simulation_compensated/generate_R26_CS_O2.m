function result = generate_R26_CS_O2(repetitions,cores)

local_job = parcluster('local');
pool = parpool(local_job, cores);

addpath('../utilities/')
addpath('../pulse_schemes/')

    
delta = -4.5e3;                     % dipolar coupling strength (linear Hz)
nu1 = 65e3;                         % rf amplitude (linear Hz)
nu_r = 10e3;                           % chemical-shift offset (linear Hz)
nucs_list = linspace(-70e3,70e3,1001); % MAS frequency list (linear Hz)  

%repetitions = 10;                   % repetitions of the pulse scheme 
scheme = scheme_R26(nu1);           
%plot_scheme(scheme);
thetam = acos(1/sqrt(3));           % magic angle (rad.)

% Fourier series coefficients parameter
npoints = 100;   % highest Fourier coefficients 
step = 0.1e-6;   % time resolution (sec.)
csflag = 1;      % include CS [0 False, 1 True]

% powder average parameter
count = 1:299;      % number of powder points
value1= 300;    
value2= 37;     
value3= 61;     
beta  = pi * count/value1;
alpha = 2*pi * rem((value2*count),value1)/value1;
gamma = 2*pi * rem((value3*count),value1)/value1;

%% definitions and allocations
tau_m = sum(scheme.tau); % modulation time (sec.)
nu_m = 1/tau_m;          % modulation frequency (linear Hz)
T = repetitions*tau_m;   % overall duration (sec.)

% definition of the spin operators 
[~,~,~, ...
 ~,~,~, ...
 ~,~,~, ...
 ~,~,~, ...
 ~,~,Iz1, ...
 ~,~,I1z, ...
 Iv] = spin_operators;

% allocation of the data arrays
signalA1 = zeros(size(nucs_list));
signalA2 = zeros(size(nucs_list));
signalB1 = zeros(size(nucs_list));
signalB2 = zeros(size(nucs_list));

% waitbar
% msg = sprintf("Up to second-order with %d powder points", count(end));
% f2 = waitbar(0,msg,'Name','Second-order',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% setappdata(f2,'canceling',0);


parfor nucs_index=1:length(nucs_list) % nur loop
%    tic;
    % waitbar
%     if getappdata(f2,'canceling')
%         break
%     end
%    waitbar(nucs_index/numel(nucs_list),f2)


    % calculate the Fourier coefficients
    [coeff,~,nu_eff,rot_ax] = sequence_get_coeff(...
        scheme.tau, ...   % duration of each pulse
        scheme.phi, ...   % rf phase of each pulse
        scheme.nu1, ...   % rf amplitude of pulse scheme (Hz linear)
        nucs_list(nucs_index), ...         % cs offset (Hz linear)
        step, ...         % time resolution (sec.)
        npoints, ...      % number of coefficients
        csflag);

    angle = 2*pi*nu_eff*T;
    axang = [rot_ax,angle];
    rotm = axang2rotm(axang);
    tilted_spins = rotm*[0,0,1]';

    h2 = calc_all_h2(nu_r,nu_m,nu_eff,T,npoints);
    %plot_coefficients(coeff)
    % datax2 is normalized: 1/sqrt(6) (2zz-xx-yy)
    [~, datax2, ~, ~, ~] = general_calc_twospin(coeff);
    % spin part: two-spin coefficients * spin operator
    Hn = zeros(4,4,5,2*npoints-1);
    for m=1:(2*npoints-1) % indexing nu_m
        for p=-2:2 % nu_eff
            for s=1:9 % spin operators
                Hn(:,:,p+3,m) = Hn(:,:,p+3,m) + datax2(m,p+3,s)*Iv(:,:,s);
            end
        end
    end
    % calculate h1
    h1 = calc_all_h1(nu_r,nu_m,nu_eff,npoints,T);

    
    for powder_index=1:length(beta) % powder loop
        scale = sin(beta(powder_index));
        % spatial part: coupling coefficients
        wn = zeros(1,5);
        for n=-2:2
            wn(n+3) = dlmn(2,n,0,thetam)*dlmn(2,0,n,beta(powder_index)) ...
                *exp(-1i*n*gamma(powder_index))*sqrt(3/2)*delta;
        end
        % full first-order effective Hamiltonian
        HamA = zeros(4,4);
        for m=1:(2*npoints-1) % indexing nu_m
            k = m-npoints;
            for n=-2:1:2
                for p=-2:2
                    %nuA = k*nu_m + n*nu_r + p*nu_eff; % linear Hz
                    HamA = HamA + Hn(:,:,p+3,m)*wn(n+3)*h1(n+3,m,p+3);
                end
            end
        end
        % second-order effective Hamiltonian
        for m1_index= 1:(2*npoints-1) % coefficients A index
            %m1 = m1_index-npoints;
            for m2_index= 1:(2*npoints-1) % coefficients B index
                %m2 = m2_index-npoints;
                for p1=-2:2
                    for p2=-2:2
                        % commutator is independent of n1 and n2
                        comm = commutator(Hn(:,:,p1+3,m1_index),Hn(:,:,p2+3,m2_index));
                        for n2=-2:1:2 % dipolar coupling B
                            for n1=-2:1:2 % dipolar coupling A
                                % resonance condition A (n1,n2)
                                % nuA = n1*nu_r+m1*nu_m+p1*nu_eff; % Hz linear
                                % resonance condition B (m1,m2)
                                % nuB = n2*nu_r+m2*nu_m+p2*nu_eff; % Hz linear
                                % Commutator of the Fourier coefficients and coupling strength
                                Ham2 = -1/2*wn(n1+3)*wn(n2+3)*comm;
                                % second-order basis function
                                %h2 = calc_h2(nuA,nuB,T);
                                Ham2 = Ham2 .* h2(n1+3,n2+3,m1_index,m2_index,p1+3,p2+3);
                                HamA = HamA + Ham2;
                            end
                        end
                    end
                end
            end
        end


        % the first and second-order Hamiltonian
        UA  = expm(-1i*2*pi*(HamA)*T);
        UAi = expm(+1i*2*pi*(HamA)*T);

        % detection
        sigmaA = I1z;
        det1   = I1z;
        det2   = Iz1;
        
         % detection in the effective frame
%           det2 = real(tilted_spins(1).*I1x ...
%                + tilted_spins(2).*I1y ...
%                + tilted_spins(3).*I1z);
%           det1 = real(tilted_spins(1).*Ix1 ...
%                + tilted_spins(2).*Iy1 ...
%                + tilted_spins(3).*Iz1);

        sigmaB = UA*sigmaA*UAi;
        signalA1(nucs_index) = signalA1(nucs_index)+trace(sigmaA*det1)*scale;
        signalA2(nucs_index) = signalA2(nucs_index)+trace(sigmaA*det2)*scale;
        signalB1(nucs_index) = signalB1(nucs_index)+trace(sigmaB*det1)*scale;
        signalB2(nucs_index) = signalB2(nucs_index)+trace(sigmaB*det2)*scale;

    end % powder loop
%toc;
end % w1 loop
%delete(f2); % remove waitbar

% normalize intensities with initial spin 1 intensity
result.signalA2 = real(signalA2/signalA1(1));
result.signalB2 = real(signalB2/signalA1(1));
result.signalB1 = real(signalB1/signalA1(1));
result.signalA1 = real(signalA1/signalA1(1));

result.nucs_list = nucs_list;
result.nu1 = nu1;
result.T = T;
end

%% functions

function h1 = calc_all_h1(nu_r,nu_m,nu_eff,npoints,T)
h1 = zeros(5,2*npoints-1,5);
for m=1:(2*npoints-1) % indexing nu_m
    k = m-npoints;
    for n=-2:1:2
        for p=-2:2
            nuA = k*nu_m + n*nu_r + p*nu_eff; % linear Hz
            h1(n+3,m,p+3) = calc_h1(nuA,T);
        end
    end
end
end

function h2 = calc_all_h2(nu_r,nu_m,nu_eff,T,npoints)
%memo_calc_h2 = memoize(@calc_h2);
h2 = zeros(5,5,2*npoints-1,2*npoints-1,5,5);
for m1_index= 1:(2*npoints-1) % coefficients A index
    m1 = m1_index-npoints;
    for m2_index= 1:(2*npoints-1) % coefficients B index
        m2 = m2_index-npoints;
        for n2=-2:1:2 % dipolar coupling B
            for n1=-2:1:2 % dipolar coupling A
                for p1=-2:2
                    for p2=-2:2
                        % resonance condition A (n1,n2)
                        nuA = n1*nu_r+p1*nu_eff+m1*nu_m; % Hz linear
                        % resonance condition B (m1,m2)
                        nuB = n2*nu_r+p2*nu_eff+m2*nu_m; % Hz linear
                        % second-order basis function
                        h2(n1+3,n2+3,m1_index,m2_index,p1+3,p2+3) = calc_h2(nuA,nuB,T);
                    end
                end
            end
        end
    end
end
end


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
% this is a dirty fix, because a single point leads to NaN and I don't know
% what is going on ...
if isnan(h)
    h = 0;
end
end
