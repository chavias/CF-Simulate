function  result = generate_C7_effective_MAS_powder_O1(repetitions)
% First and second-order effective Hamiltonian for postC7 recoupling
% as a function of the MAS frequency nur for fixed rf amplitude nu1

% WITHOUT CHEMICAL SHIFT

addpath('../../utilities/', "../../pulse_schemes/");
% experimental parameter
delta = -4.5e3;                      % dipolar coupling strength (linear Hz)
nu1 = 70e3;                         % rf amplitude (linear Hz)
nucs = 0;                           % chemical-shift offset (linear Hz)
nur_list = linspace(2e3,21e3,1001); % MAS frequency list (linear Hz)  

scheme = scheme_C7(nu1);       
%plot_scheme(scheme);               
thetam = acos(1/sqrt(3));           % magic angle (rad.)

% Fourier series coefficients parameter
npoints = 100;   % highest Fourier coefficients 
step = 0.1e-6;   % time resolution (sec.)
csflag = 0;      % include CS [0 False, 1 True]

% powder average parameter
count = 1:299; 
value1= 300;
%value2= 37;
value3= 61;
beta  = pi * count/value1;
%alpha = 2*pi * rem((value2*count),value1)/value1;
gamma = 2*pi * rem((value3*count),value1)/value1;

% definitions and allocations
tau_m = sum(scheme.tau); % modulation time (sec.)
nu_m = 1/tau_m;          % modulation frequency (linear Hz)
T = repetitions*tau_m;   % overall duration (sec.)

% numerical data
% dataA = squeeze(data(:,1,index_T(T,time_resolution)))./data(1,1,1);
% dataB = squeeze(data(:,2,index_T(T,time_resolution)))./data(1,1,1);

% definition of the spin operators 
[~,~,~, ...
 ~,~,~, ...
 ~,~,~, ...
 ~,~,~, ...
 ~,~,Iz1, ...
 ~,~,I1z,Iv] = spin_operators;


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

% Second-order terms

% data allocation
signal2A1=zeros(size(nur_list));
signal2A2=zeros(size(nur_list));
signal2B1=zeros(size(nur_list));
signal2B2=zeros(size(nur_list));


%% idea precalculate the coupling and store it in an array

% calculate spatial part 
[wn,scale] = calc_coupling(beta,gamma,delta,thetam);

% calculat spin part
Hn = calc_spin_part(datax2,Iv,npoints);

% calculate all h1
h1 = calc_all_h1(nu_m,npoints,nur_list,T);

parfor nur_index=1:length(nur_list) % nur loop
    %tic;

    for powder_index=1:length(beta) % powder loop
        
        % possible optimization: could be outside the loop
        % spin part: two-spin coefficients * spin operator
        
        % full Hamiltonian
        HamA = zeros(4,4);
        for m=1:(2*npoints-1) % indexing wm
            for n=-2:1:2
                % effective frequency is not present
                %k = m-npoints; 
                % nuA = k*nu_m + n*nur_list(nur_index); % linear in Hz
                HamA = HamA + Hn(:,:,m).*wn(n+3,powder_index).*h1(m,n+3,nur_index);
            end
        end


    % the first and second-order Hamiltonian
    UA  = expm(-1i*2*pi*(HamA)*T);
    UAi = expm(+1i*2*pi*(HamA)*T);

    % detection
    sigmaA = I1z;
    det1   = I1z;
    det2   = Iz1;
    sigmaB = UA*sigmaA*UAi;

    signal2A1(nur_index) = signal2A1(nur_index)+trace(sigmaA*det1)*scale(powder_index);
    signal2A2(nur_index) = signal2A2(nur_index)+trace(sigmaA*det2)*scale(powder_index);
    signal2B1(nur_index) = signal2B1(nur_index)+trace(sigmaB*det1)*scale(powder_index);
    signal2B2(nur_index) = signal2B2(nur_index)+trace(sigmaB*det2)*scale(powder_index);
    end  % powder loop
%toc;
end % w1 loop
%delete(f2); % remove waitbar

% normalize intensities with initial spin 1 intenstiy
result.signalA2 = real(signal2A2/signal2A1(1));
result.signalB2 = real(signal2B2/signal2A1(1));
result.signalB1 = real(signal2B1/signal2A1(1));
result.signalA1 = real(signal2A1/signal2A1(1));

result.nur_list = nur_list;
result.nu1 = nu1;
result.T = T;
end 


%% functions

% 
% function h2 = calc_all_h2(nur_list, nu_m, T, npoints)
% 
% %memo_calc_h2 = memoize(@calc_h2);
% 
% h2 = zeros(5, 5, 2*npoints-1, 2*npoints-1);
% for nur_index = 1:length(nur_list) % nur loop
%     tic;
%     for m1_index = 1:(2*npoints-1) % coefficients A index
%         m1 = m1_index - npoints;
%         for m2_index = 1:(2*npoints-1) % coefficients B index
%             m2 = m2_index - npoints;
%             for n2 = -2:1:2 % dipolar coupling B
%                 for n1 = -2:1:2 % dipolar coupling A
%                     % Compute unique resonance conditions only
%                     unique_resonances = unique(bsxfun(@plus, n1 * nur_list(nur_index), m1 * nu_m) + ...
%                                                bsxfun(@plus, n2 * nur_list(nur_index), m2 * nu_m)');
%                     for resonance_index = 1:length(unique_resonances)
%                         nuA = unique_resonances(resonance_index);
%                         % Compute resonance condition B only once
%                         nuB = n2 * nur_list(nur_index) + m2 * nu_m; % Hz linear
%                         % Check if nuB was computed before
%                         if nuB >= nuA
%                             % Compute second-order basis function
%                             h2(n1+3, n2+3, m1_index, m2_index, nur_index) = calc_h2(nuA, nuB, T);
%                         else
%                             % Use the symmetry property
%                             h2(n1+3, n2+3, m1_index, m2_index, nur_index) = h2(n2+3, n1+3, m2_index, m1_index, nur_index);
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     toc;
% end
% end

% function h2 = calc_all_h2(nur_list, nu_m, T, npoints)
%     % Initialize h2 cell array
%     h2_cell = cell(numel(nur_list), 1);
% 
%     % Loop over nur_list
%     parfor nur_index = 1:numel(nur_list)
%         % Initialize h2 matrix for this iteration
%         h2_iter = zeros(5, 5, 2*npoints-1, 2*npoints-1);
% 
%         % Loop over coefficients A index
%         for m1_index = 1:(2*npoints-1)
%             m1 = m1_index - npoints;
% 
%             % Loop over coefficients B index
%             for m2_index = 1:(2*npoints-1)
%                 m2 = m2_index - npoints;
% 
%                 % Loop over dipolar coupling B
%                 for n2 = -2:2
%                     % Loop over dipolar coupling A
%                     for n1 = -2:2
%                         % Resonance condition A (n1, n2)
%                         nuA = n1 * nur_list(nur_index) + m1 * nu_m;
% 
%                         % Resonance condition B (m1, m2)
%                         nuB = n2 * nur_list(nur_index) + m2 * nu_m;
% 
%                         % Calculate second-order basis function
%                         h2_iter(n1+3, n2+3, m1_index, m2_index) = calc_h2(nuA, nuB, T);
%                     end
%                 end
%             end
%         end
% 
%         % Store the calculated h2 matrix for this iteration
%         h2_cell{nur_index} = h2_iter;
%     end
% 
%     % Convert cell array to h2 array
%     h2 = cat(5, h2_cell{:});
% end


function h2 = calc_all_h2(nur_index,nu_m,T,npoints)
h2 = zeros(5,5,2*npoints-1,2*npoints-1);
    for m1_index= 1:(2*npoints-1) % coefficients A index
        m1 = m1_index-npoints;
        for m2_index= 1:(2*npoints-1) % coefficients B index
            m2 = m2_index-npoints;
            for n2=-2:1:2 % dipolar coupling B
                for n1=-2:1:2 % dipolar coupling A
                    % resonance condition A (n1,n2)
                    nuA = n1*nur_list(nur_index)+m1*nu_m; % Hz linear
                    % resonance condition B (m1,m2)
                    nuB = n2*nur_list(nur_index)+m2*nu_m; % Hz linear
                    % second-order basis function
                    h2(n1+3,n2+3,m1_index,m2_index,nur_index) = calc_h2(nuA,nuB,T);
                end
            end
        end
    end
end

function comm = calc_all_comm(Hn,npoints)
comm = zeros(4,4,2*npoints-1,2*npoints-1);
for m1_index= 1:(2*npoints-1) % coefficients A index
    for m2_index= 1:(2*npoints-1) % coefficients B index
        comm(:,:,m1_index,m2_index) = commutator(Hn(:,:,m1_index),Hn(:,:,m2_index));
    end
end
end


function [wn,scale] = calc_coupling(beta,gamma,delta,thetam)
% spatial part: coupling coefficients
wn = zeros(1,5,numel(beta));
scale = zeros(1,numel(beta));
for powder_index=1:length(beta) % powder loop
    scale(powder_index)=sin(beta(powder_index));
    for n=-2:2
        wn(n+3,powder_index) = dlmn(2,n,0,thetam)*dlmn(2,0,n,beta(powder_index)) ...
            *exp(-1i*n*gamma(powder_index))*sqrt(3/2)*delta;
    end
end
end

function Hn = calc_spin_part(datax2,Iv,npoints)
Hn = zeros(4,4,2*npoints-1);
for m=1:(2*npoints-1) % indexing wm
    for s=1:9 % spin operators
        Hn(:,:,m) = Hn(:,:,m) + datax2(m,3,s)*Iv(:,:,s);
    end
end
end

function h1 = calc_all_h1(nu_m,npoints,nur_list,T)
h1 = zeros(2*npoints-1,5,numel(nur_list));
for nur_index=1:length(nur_list) % nur loop
        for m=1:(2*npoints-1) % indexing wm
            for n=-2:1:2
                % effective frequency is not present
                k = m-npoints; 
                nuA = k*nu_m + n*nur_list(nur_index); % linear in Hz
                h1(m,n+3,nur_index) =  calc_h1(nuA,T);
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

