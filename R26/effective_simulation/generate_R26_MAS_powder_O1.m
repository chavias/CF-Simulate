function result = generate_R26_MAS_powder_O1(repetitions)

delta = -4.5e3;                     % dipolar coupling strength (linear Hz)
nu1 = 65e3;                         % rf amplitude (linear Hz)
nucs = 0;                           % chemical-shift offset (linear Hz)
nur_list = linspace(2e3,21e3,1001); % MAS frequency list (linear Hz)  

%repetitions = 10;                   % repetitions of the pulse scheme 
scheme = scheme_R26(nu1);           
%plot_scheme(scheme);
thetam = acos(1/sqrt(3));           % magic angle (rad.)

% Fourier series coefficients parameter
npoints = 120;   % highest Fourier coefficients 
step = 0.1e-6;   % time resolution (sec.)
csflag = 1;      % include CS [0 False, 1 True]

% powder average parameter
count = 1:1;      % number of powder points
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
% index_T(T,time_resolution)
% % numerical data
% dataA = squeeze(data(:,1,index_T(T,time_resolution)))./data(1,1,1);
% dataB = squeeze(data(:,2,index_T(T,time_resolution)))./data(1,1,1);
% 
% plot(nur_list/1e3,dataB)

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
[coeff,~,nu_eff,rot_ax] = sequence_get_coeff(...
    scheme.tau, ...   % duration of each pulse
    scheme.phi, ...   % rf phase of each pulse
    scheme.nu1, ...   % rf amplitude of pulse scheme (Hz linear)
    nucs, ...         % cs offset (Hz linear)
    step, ...         % time resolution (sec.)
    npoints, ...      % number of coefficients
    csflag);

angle = 2*pi*nu_eff*T;
axang = [rot_ax,angle];
rotm = axang2rotm(axang);
tilted_spins = rotm*[0,0,1]';

%plot_coefficients(coeff)
% datax2 is normalized: 1/sqrt(6) (2zz-xx-yy)
[~, datax2, ~, ~, ~] = general_calc_twospin(coeff);

% data allocation
signal2A1=zeros(size(nur_list));
signal2A2=zeros(size(nur_list));
signal2B1=zeros(size(nur_list));
signal2B2=zeros(size(nur_list));

% calculate spatial part 
[wn,scale] = calc_coupling(beta,gamma,delta,thetam);

% calculat spin part
Hn = calc_spin_part(datax2,Iv,npoints);


for nur_index=1:length(nur_list) % nur loop

    % calculate all h1
    h1 = calc_all_h1(nu_m,npoints,nur_list,nur_index,T,nu_eff);

    for powder_index=1:length(beta) % powder loop
        
        % possible optimization: could be outside the loop
        % spin part: two-spin coefficients * spin operator

        % full Hamiltonian
        HamA = zeros(4,4);
        for m=1:(2*npoints-1) % indexing wm
            for p=-2:2
            for n=-2:1:2
                % effective frequency is not present
                %k = m-npoints; 
                % nuA = k*nu_m + n*nur_list(nur_index); % linear in Hz
                HamA = HamA + Hn(:,:,m,p+3).*wn(n+3,powder_index).*h1(m,n+3,p+3);
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
    sigmaB = UA*sigmaA*UAi;

    signal2A1(nur_index) = signal2A1(nur_index)+trace(sigmaA*det1)*scale(powder_index);
    signal2A2(nur_index) = signal2A2(nur_index)+trace(sigmaA*det2)*scale(powder_index);
    signal2B1(nur_index) = signal2B1(nur_index)+trace(sigmaB*det1)*scale(powder_index);
    signal2B2(nur_index) = signal2B2(nur_index)+trace(sigmaB*det2)*scale(powder_index);
    end  % powder loop
end % w1 loop


% normalize intensities with initial spin 1 intensity
result.signalA2 = real(signal2A2/signal2A1(1));
result.signalB2 = real(signal2B2/signal2A1(1));
result.signalB1 = real(signal2B1/signal2A1(1));
result.signalA1 = real(signal2A1/signal2A1(1));

result.nur_list = nur_list;
result.nu1 = nu1;
result.T = T;
end


function comm = calc_all_comm(Hn,npoints)
comm = zeros(4,4,2*npoints-1,2*npoints-1,5,5);
for m1_index= 1:(2*npoints-1) % coefficients A index
    for m2_index= 1:(2*npoints-1) % coefficients B index
        for p1=-2:2
            for p2=-2:2
                comm(:,:,m1_index,m2_index,p1+3,p2+3) = commutator(Hn(:,:,m1_index,p1+3),Hn(:,:,m2_index,p2+3));
            end
        end
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
Hn = zeros(4,4,2*npoints-1,5);
for p = -2:2
for m=1:(2*npoints-1) % indexing wm
    for s=1:9 % spin operators
        Hn(:,:,m,p+3) = Hn(:,:,m,p+3) + datax2(m,p+3,s)*Iv(:,:,s);
    end
end
end
end

function h1 = calc_all_h1(nu_m,npoints,nur_list,nur_index,T,nu_eff)
h1 = zeros(2*npoints-1,5,5);
    for p=-2:2
        for m=1:(2*npoints-1) % indexing wm
            for n=-2:1:2
                % effective frequency is not present
                k = m-npoints; 
                nuA = k*nu_m + n*nur_list(nur_index)+p*nu_eff; % linear in Hz
                h1(m,n+3,p+3) =  calc_h1(nuA,T);
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