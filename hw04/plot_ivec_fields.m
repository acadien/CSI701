%close all;
clear all;
load ixvec_fld.dat;
load ixpts.dat;
load iyvec_fld.dat;
load iypts.dat;

X=ixpts;
Y=iypts;
U=ixvec_fld;
V=iyvec_fld;
vfld_mag=(U.^2+V.^2).^(1/2);
Q=100;

figure;
contour(X,Y,vfld_mag,20);
colorbar;
title('||Gradient Vector Field Contour|| - Finite Element Method');

figure;
mesh(X,Y,vfld_mag);
colorbar;
title('||Gradient Vector Field Surface|| - Finite Element Method');

figure;
surf(X,Y,vfld_mag);
colorbar;
title('||Gradient Vector Field Surface|| - Finite Element Method');

