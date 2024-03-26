function ret = dlmn(l, m, n, beta)

% calculate wigner rotation matrix elements dlmn for
% a given angle beta.
%
% Parameters:
%	l	rank of tensor
%	m,n	matrixelements
%	beta	angle in rad !!!!!!!!!!!!
%

ret  = 0;
for v=max(0,n-m):min(l-m,l+n)
    ret = ret + (-1) ^ (v+m-n) * ...
	(cos(beta/2).^(2*l+n-m-2*v) .* (sin(beta/2)).^(2*v+m-n) / ...
	(fac(l+n-v) * fac(l-m-v) * fac(v) * fac(v+m-n)));
end
ret = sqrt(fac(l+m)*fac(l-m)*fac(l+n)*fac(l-n)) * ret;

