%close all;
load fld.dat;
load pts.dat;
load vec_field.dat;

scl_fld=fld;
X=pts(:,1);
Y=pts(:,2);
U=vec_field(:,1);
V=vec_field(:,2);

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
            %vfld_mag(i,j)=0;
            %sfld(i,j)=0;
            %sfld_mag(i,j)=0;
        end
    end
end

figure;
contour(xi,yi,vfld_mag,20);
colorbar;
title('||Gradient Vector Field||');

figure;
contour(xi,yi,sfld,20);
colorbar;
title('Scalar Field Contour');

%figure;
%contour(xi,yi,sfld_mag,20);
%colorbar;
%title('Scalar Field Contour of Magnitude abs(scalar-field)');

