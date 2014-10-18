close all;
%clear all;
if(1)
load ptcls_x.dat;
load ptcls_y.dat;
load ptcls_vx.dat;
load ptcls_vy.dat;
end

prefix='./pics/ptcl_stage';
suffix='.eps';

x=0.3*cos(linspace(0,2*pi,100));
y=0.3*sin(linspace(0,2*pi,100));


fig=figure('visible','off');
for q=1:(squeeze(size(ptcls_x(:,1))))
    %pause(0.1)
    a=sqrt(ptcls_vx(q,:).^2+ptcls_vy(q,:).^2);
    frm=scatter(ptcls_x(q,:),ptcls_y(q,:),30*a+0.1,a,'fillon');
    hold on;
    plot(x,y,'k');
    %scatter(ptcls_x(q,:),ptcls_y(q,:),,'k');
    hold off;
    caxis('manual');
    axis([-1,1,-1,1]);
    colorbar;
    num=sprintf('%3.3d',q);
    flnm=strcat(prefix,num);
    saveas(frm,flnm,'jpeg');
end
