close all;
load act_fld.dat;
load ixpts.dat;
load iypts.dat;

X=ixpts;
Y=iypts;
Q=400;
Z=zeros(Q,100,100);
for k=1:Q
    Z(k,1:100,1:100)=act_fld(((k-1)*100+1):(k*100),1:100);
end

fig=figure('visible','off');

prefix='./pics/shade_stage';
for i=1:400
    frm = surf(X,Y,squeeze(Z(i,:,:)));
    axis([-1 1 -1 1 0 1]); 
    view([0,0,1]);
    shading interp;
    num=sprintf('%3.3d',i);
    flnm=strcat(prefix,num);
    saveas(frm,flnm,'jpeg');
end
