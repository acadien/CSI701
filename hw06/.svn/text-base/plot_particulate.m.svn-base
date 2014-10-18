close all;
if(1)
load ptcls_x.dat;
load ptcls_y.dat;
load ptcls_vx.dat;
load ptcls_vy.dat;
end

a1=2; a2=50; a3=100; a4=160;

X1=ptcls_x(a1,:);
X2=ptcls_x(a2,:);
X3=ptcls_x(a3,:);
X4=ptcls_x(a4,:);
Y1=ptcls_y(a1,:);
Y2=ptcls_y(a2,:);
Y3=ptcls_y(a3,:);
Y4=ptcls_y(a4,:);

ZX1=ptcls_vx(a1,:);
ZX2=ptcls_vx(a2,:);
ZX3=ptcls_vx(a3,:);
ZX4=ptcls_vx(a4,:);
ZY1=ptcls_vy(a1,:);
ZY2=ptcls_vy(a2,:);
ZY3=ptcls_vy(a3,:);
ZY4=ptcls_vy(a4,:);


figure();

sp(1)=subplot(2,2,1);
scatter(X1,Y1,10,sqrt(ZX1.^2+ZY1.^2),'fillon');
caxis('manual');
axis([-1,1,-1,1,0,1]);
title('Initial Field, Timestep 2, color=|v|');
view([0,0,1]);

sp(2)=subplot(2,2,2);
scatter(X2,Y2,10,sqrt(ZX2.^2+ZY2.^2),'fillon');
caxis('manual');
axis([-1,1,-1,1,0,1]);
title('Field at Timestep 50, color=|v|');
view([0,0,1]);

sp(3)=subplot(2,2,3);
scatter(X3,Y3,10,sqrt(ZX3.^2+ZY3.^2),'fillon');
caxis('manual');
axis([-1,1,-1,1,0,1]);
title('Field at Timestep 100, color=|v|');
view([0,0,1]);

sp(4)=subplot(2,2,4);
scatter(X4,Y4,10,sqrt(ZX4.^2+ZY4.^2),'fillon');
caxis('manual');
axis([-1,1,-1,1,0,1]);
title('Field at Timestep 160, color=|v|');
view([0,0,1]);

h=colorbar;
set(h,'Position',[0.9314 0.1 0.0281 0.8150])
