function [t1, data, nu_eff, rot_ax] = sequence_iframe9_time(tau,phi,nu1,nucs,step,csflag)
% calculates the interaction frame of a general sequence
% with a pulse time vector tau and a phase vector phi.
% tau:   list of pulses lengths (sec.)
% phi:   list of phase values (degree.)
% nu1:   rf field amplitude (linear Hz)
% nucs   chemical-shift offset (linear Hz)
% step:  time resolution for calculation (sec.)
% csflag: include chemical shift during pulse

if nargin < 6
  csflag = 1;
end

w1=nu1*2*pi;   % angular Hz
wcs=nucs*2*pi; % angular Hz

t1=0:step:(sum(tau)+step/2);
ampl1 = ones(size(t1));

for k=length(tau):-1:1
  ampl1(t1<(sum(tau(1:k))-step/2))=exp(1i*phi(k)/180*pi);
end
ampl1(end)=ampl1(1);

%subplot(1,2,1)
%plot(t1,real(ampl1),'xb',t1,imag(ampl1),'rx');
%subplot(1,2,2)
%plot(t1,abs(ampl1),'xb',t1,angle(ampl1)*180/pi,'rx');

ampl = w1*ampl1;

Sx = 1/2*[0 1;1 0];
Sy = 1i/2*[0 -1;1 0];
Sz = 1/2*[1 0;0 -1];
Se = [1 0;0 1];

% for performance it would be better to implement the solution for the
% coefficients directly
U = zeros(2,2,length(t1));
for k1=1:length(t1)
  if csflag
    U(:,:,k1) = expm(1i*(wcs*Sz+real(ampl(k1))*Sx+imag(ampl(k1))*Sy)*step);
  else
    U(:,:,k1) = expm(1i*(real(ampl(k1))*Sx+imag(ampl(k1))*Sy)*step);
  end
end
xx = zeros(size(t1));
xy = zeros(size(t1));
xz = zeros(size(t1));
yx = zeros(size(t1));
yy = zeros(size(t1));
yz = zeros(size(t1));
zx = zeros(size(t1));
zy = zeros(size(t1));
zz = zeros(size(t1));
Up=Se;
xx(1)=1;
yy(1)=1;
zz(1)=1;
for k1=2:length(t1)
  Up=Up*U(:,:,k1-1);
  sigmax=Up*Sx*Up';
  sigmay=Up*Sy*Up';
  sigmaz=Up*Sz*Up';
  xx(k1) = 2*trace(sigmax*Sx);
  xy(k1) = 2*trace(sigmax*Sy);
  xz(k1) = 2*trace(sigmax*Sz);
  yx(k1) = 2*trace(sigmay*Sx);
  yy(k1) = 2*trace(sigmay*Sy);
  yz(k1) = 2*trace(sigmay*Sz);
  zx(k1) = 2*trace(sigmaz*Sx);
  zy(k1) = 2*trace(sigmaz*Sy);
  zz(k1) = 2*trace(sigmaz*Sz);
end

data=zeros(9,length(t1));
data(1,:) = xx;
data(2,:) = xy;
data(3,:) = xz;
data(4,:) = yx;
data(5,:) = yy;
data(6,:) = yz;
data(7,:) = zx;
data(8,:) = zy;
data(9,:) = zz;

R=zeros(3,3);
R(1,1)=data(1,end);
R(1,2)=data(2,end);
R(1,3)=data(3,end);
R(2,1)=data(4,end);
R(2,2)=data(5,end);
R(2,3)=data(6,end);
R(3,1)=data(7,end);
R(3,2)=data(8,end);
R(3,3)=data(9,end);
R=real(R);

r=real(rotm2axang(R)); 
rot_ax=r(1:3);
theta=r(4);
nu_eff=theta/(2*pi*sum(tau));
