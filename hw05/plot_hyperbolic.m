
load act_fld.dat;
load ixpts.dat;
load iypts.dat;

X=ixpts;
Y=iypts;
Z1(1:100,1:100)=act_fld(((2-1)*100+1):(2*100),1:100);
Z2(1:100,1:100)=act_fld(((100-1)*100+1):(100*100),1:100);
Z3(1:100,1:100)=act_fld(((200-1)*100+1):(200*100),1:100);
Z4(1:100,1:100)=act_fld(((300-1)*100+1):(300*100),1:100);

figure;

sp(1)=subplot(2,2,1);
surf(X,Y,Z1);
caxis('manual');
axis([-1,1,-1,1,0,1]);
title('Initial Field, Timestep 20/8000');
view([0,0,1]);

sp(2)=subplot(2,2,2);
surf(X,Y,Z2);
caxis('manual');
axis([-1,1,-1,1,0,1]);
title('Field at Timestep 2000/8000');
view([0,0,1]);

sp(3)=subplot(2,2,3);
surf(X,Y,Z3);
caxis('manual');
axis([-1,1,-1,1,0,1]);
title('Field at Timestep 4000/8000');
view([0,0,1]);

sp(4)=subplot(2,2,4);
surf(X,Y,Z4);
caxis('manual');
axis([-1,1,-1,1,0,1]);
title('Field at Timestep 6000/8000');
view([0,0,1]);

h=colorbar;
set(h,'Position',[0.9314 0.1 0.0281 0.8150])
