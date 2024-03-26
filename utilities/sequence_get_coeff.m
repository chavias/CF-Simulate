function [coeff, nu_m, nu_eff, rot_ax] = sequence_get_coeff(tau,phi,nu1,nucs,step,npoints,csflag)
% calculate the Fourier coefficients of the interaction frame of CPMG
%
% Input
%   tau:   length of pulses
%   phi:   phase values
%   nu1:   rf field amplitude (linear Hz)
%   nucs:  chemical shift (linear Hz)
%   step:  time resolution for calculation
%   npoints: highest Fourier coefficients
%   csflag: include chemical shift during pulse
%
% The highest Fourier coefficient in nu_m is npoints-1!
%

% round tau values such that they are multiples of step;
tau = round(tau/step)*step;

[time, data, nu_eff, rot_ax] = sequence_iframe9_time(tau,phi,nu1,nucs,step,csflag);
nu_m = 1/(sum(tau));
nu_eff = min(nu_eff,nu_m-nu_eff);
% phi = atan2(rot_ax(2),rot_ax(1))*180/pi;

%fprintf('nu_{mod} = %8.3f kHz\n',nu_m/(1000));
%fprintf('nu_{res} = %6.1f Hz\n',1/(sum(tau)));
%fprintf('nu_{eff} = %8.3f Hz\n',nu_eff);
%fprintf('theta = %8.1f degree\n',theta);
%fprintf('phi = %8.1f degree\n',phi);

% special treatment for effective field
if (abs(nu_eff) > 0 && abs(nu_eff-nu_m) > 0)

    % algin with effective field axis
    theta = acos(rot_ax(3));            % y-rotation (rad.)
    psi   = atan2(rot_ax(2),rot_ax(1)); % z-rotation (rad.)
    data1=zeros(size(data));
    data1(1:3,:) = Ry(Rz(data(1:3,:),-psi),-theta);
    data1(4:6,:) = Ry(Rz(data(4:6,:),-psi),-theta);
    data1(7:9,:) = Ry(Rz(data(7:9,:),-psi),-theta);
    
    % rotate with effective field 
    data2=zeros(size(data));
    data2([3 6 9],:)=data1([3 6 9],:);
    co=cos(2*pi*nu_eff*time);
    si=sin(2*pi*nu_eff*time);
    data2(1,:) = data1(1,:).*co-data1(2,:).*si;
    data2(2,:) = data1(2,:).*co+data1(1,:).*si;
    data2(4,:) = data1(4,:).*co-data1(5,:).*si;
    data2(5,:) = data1(5,:).*co+data1(4,:).*si;
    data2(7,:) = data1(7,:).*co-data1(8,:).*si;
    data2(8,:) = data1(8,:).*co+data1(7,:).*si;
    
    % fourier transform
    scale=length(data2)-1;
    dataFT2=zeros(9,scale);
    for k=1:9
        dataFT2(k,:) = fftshift(fft(data2(k,1:end-1)))/scale;
    end
    
    % rotate with effective field back
    dataFT1 = zeros(scale,3,9);
    dataFT1(:,2,3) = dataFT2(3,:);
    dataFT1(:,2,6) = dataFT2(6,:);
    dataFT1(:,2,9) = dataFT2(9,:);

    dataFT1(:,1,1) = 0.5*(dataFT2(1,:)+1i*dataFT2(2,:));
    dataFT1(:,3,1) = 0.5*(dataFT2(1,:)-1i*dataFT2(2,:));
    dataFT1(:,1,2) = 0.5*(dataFT2(2,:)-1i*dataFT2(1,:));
    dataFT1(:,3,2) = 0.5*(dataFT2(2,:)+1i*dataFT2(1,:));

    dataFT1(:,1,4) = 0.5*(dataFT2(4,:)+1i*dataFT2(5,:));
    dataFT1(:,3,4) = 0.5*(dataFT2(4,:)-1i*dataFT2(5,:));
    dataFT1(:,1,5) = 0.5*(dataFT2(5,:)-1i*dataFT2(4,:));
    dataFT1(:,3,5) = 0.5*(dataFT2(5,:)+1i*dataFT2(4,:));

    dataFT1(:,1,7) = 0.5*(dataFT2(7,:)+1i*dataFT2(8,:));
    dataFT1(:,3,7) = 0.5*(dataFT2(7,:)-1i*dataFT2(8,:));
    dataFT1(:,1,8) = 0.5*(dataFT2(8,:)-1i*dataFT2(7,:));
    dataFT1(:,3,8) = 0.5*(dataFT2(8,:)+1i*dataFT2(7,:));
    
    % align with with the original axis
    coeff=zeros(size(dataFT1));
    for k=1:3
        coeff(:,k,1:3) = Rz(Ry(squeeze(dataFT1(:,k,1:3)),theta),psi);
        coeff(:,k,4:6) = Rz(Ry(squeeze(dataFT1(:,k,4:6)),theta),psi);
        coeff(:,k,7:9) = Rz(Ry(squeeze(dataFT1(:,k,7:9)),theta),psi);
    end
else
    scale=length(data)-1;
    coeff = zeros(scale,3,9);
    for k=1:9
        coeff(:,2,k) = fftshift(fft(data(k,1:end-1)))/scale;
    end
end

zz = (-scale/2):(scale/2-1);
z1 = find(zz==-(npoints-1));
z2 = find(zz==(npoints-1));
coeff = coeff(z1:z2,:,:);
  