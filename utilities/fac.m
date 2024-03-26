function ret = fac(m)

% calculates the factorial of m if
% m is a positive integer number.
% Otherwise it returns one

if(m<1 | m ~= round(m))
  ret=1;
else
  ret=1;
  for z=m:-1:1
    ret = ret * z;
  end
end
