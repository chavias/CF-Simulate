function [datax1, datax2, datas0, datas1, datas2] = general_calc_twospin(coeff)
% calculate two spin coefficients

datax1 = coeff(:,:,7:9);
datax2 = zeros(2*length(datax1)-1,5,9);
% scale=prod(size(datax1(:,:,1)));

for k=1:3
  for j=1:3
    datax2(:,:,(k-1)*3+j)=conv2(datax1(:,:,k),datax1(:,:,j));
  end
end
fmin = floor(length(datax1)/2);
datax2=datax2(fmin+1:fmin+length(datax1),:,:);

datax2(fmin+1,3,1)=datax2(fmin+1,3,1)-1/3;
datax2(fmin+1,3,5)=datax2(fmin+1,3,5)-1/3;
datax2(fmin+1,3,9)=datax2(fmin+1,3,9)-1/3;
datax2=datax2*3/sqrt(6); % normalized

% spherical coefficients
datas0=zeros(length(datax2),5,1);
datas1=zeros(length(datax2),5,3);
datas2=zeros(length(datax2),5,5);


datas0(:,:,1) = -1/sqrt(3)*(datax2(:,:,1)+datax2(:,:,5)+datax2(:,:,9));

datas1(:,:,1) = -1/2*(datax2(:,:,7)-datax2(:,:,3)-1i*(datax2(:,:,8)-datax2(:,:,6)));
datas1(:,:,2) = -1i/sqrt(2)*(datax2(:,:,2)-datax2(:,:,4));
datas1(:,:,3) = -1/2*(datax2(:,:,7)-datax2(:,:,3)+1i*(datax2(:,:,8)-datax2(:,:,6)));

datas2(:,:,1) = +(1/2*datax2(:,:,1)-1/2*datax2(:,:,5)-1i/2*datax2(:,:,2)-1i/2*datax2(:,:,4));
datas2(:,:,2) = +(1/2*datax2(:,:,7)-1i/2*datax2(:,:,8)+1/2*datax2(:,:,3)-1i/2*datax2(:,:,6));
datas2(:,:,3) = 1/sqrt(6)*(2*datax2(:,:,9)-datax2(:,:,1)-datax2(:,:,5));
datas2(:,:,4) = -(1/2*datax2(:,:,7)+1i/2*datax2(:,:,8)+1/2*datax2(:,:,3)+1i/2*datax2(:,:,6));
datas2(:,:,5) = +(1/2*datax2(:,:,1)-1/2*datax2(:,:,5)+1i/2*datax2(:,:,2)+1i/2*datax2(:,:,4));

