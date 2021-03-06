%close all;
clear all;
load scl_fld.dat;
load vec_fld.dat;
load pts.dat;

X=pts(:,1);
Y=pts(:,2);
U=vec_fld(:,1);
V=vec_fld(:,2);
%quiver(X,Y,U,V,'filled','LineWidth',2)
%axis([-1,1,-1,1])

a=[X,Y,zeros(1413,1)];
b=[X+U,Y+V,zeros(1413,1)];
nm=cross(a,b);
nm=nm(:,3);
mag=sqrt(U.^2+V.^2);
Q=500;
ti=linspace(-1,1,Q);
[xi,yi]=meshgrid(ti,ti);

vfld_mag=griddata(X,Y,mag,xi,yi);
%vfld_nrm=griddata(X,Y,nm,xi,yi);
sfld=griddata(X,Y,scl_fld,xi,yi);
sfld_mag=griddata(X,Y,abs(scl_fld),xi,yi);

for(i=1:Q)
    for(j=1:Q)
        if(xi(i,j)^2+yi(i,j)^2<0.3^2)
            vfld_mag(i,j)=0;
            sfld(i,j)=0;
            %sfld_mag(i,j)=0;
        end
    end
end
%plot(1:1413,

close all;
figure;
contour(xi,yi,vfld_mag,20);
colorbar;
title('||Gradient Vector Field|| - Elliptical Solution');

figure;
contour(xi,yi,sfld,20);
colorbar;
title('Scalar Field Contour - Elliptical Solution');

