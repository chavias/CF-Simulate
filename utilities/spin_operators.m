function [Ix,Iy,Iz, ...
          Ixx,Ixy,Ixz, ...
          Iyx,Iyy,Iyz, ...
          Izx,Izy,Izz, ...
          Ix1,Iy1,Iz1, ...
          I1x,I1y,I1z, ...
          Iv] = spin_operators
% defines all one and two spin operator in xyz basis

Ix = 1/2 *[0  1;1  0];
Iy = 1i/2*[0 -1;1  0];
Iz = 1/2 *[1  0;0 -1];

Ixx = kron(Ix,Ix);
Ixy = kron(Ix,Iy);
Ixz = kron(Ix,Iz);

Iyx = kron(Iy,Ix);
Iyy = kron(Iy,Iy);
Iyz = kron(Iy,Iz);

Izx = kron(Iz,Ix);
Izy = kron(Iz,Iy);
Izz = kron(Iz,Iz);

Iv = zeros(4,4,9);
Iv(:,:,1) = Ixx;
Iv(:,:,2) = Ixy;
Iv(:,:,3) = Ixz;

Iv(:,:,4) = Iyx;
Iv(:,:,5) = Iyy;
Iv(:,:,6) = Iyz;

Iv(:,:,7) = Izx;
Iv(:,:,8) = Izy;
Iv(:,:,9) = Izz;

Ix1 = kron(Ix,eye(2));
Iy1 = kron(Iy,eye(2));
Iz1 = kron(Iz,eye(2));

I1x = kron(eye(2),Ix);
I1y = kron(eye(2),Iy);
I1z = kron(eye(2),Iz);
end