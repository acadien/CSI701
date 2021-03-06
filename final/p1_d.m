close all;
clear;
L=10; N=15;
xx=linspace(0,L,N);
u_f=zeros(N,1);
u_f(1)=1.0;
u_f(N)=2.0;
u_b=u_f;

%Equation Constants
k=2;
delx2=(xx(2)-xx(1))^2;
delt=delx2/3; %Set via stability requirements for Explicit scheme.
g=k*delt/delx2;

%A Matrix for progressing the implicit schema
A=eye(N)*(1+2*g);
a=eye(N-1)*(-g);
A(2:N,1:N-1)=A(2:N,1:N-1)+a;
A(1:N-1,2:N)=A(1:N-1,2:N)+a;
A(1,1)=1; %Boundary Conditions
A(1,2)=0;
A(N,N)=1;
A(N,N-1)=0;

k=0;
for t=1:100
    
    %Explicit Scheme
    for i=1:N
        if i~=1 && i~=N
            u_f(i)=g*(u_f(i+1)-2*u_f(i)+u_f(i-1))+u_f(i);           
        end
    end
    
    %Implicit Scheme
    u_b=A\u_b;
    
    if t==2 || t==10 || t==25 || t==100
        k=k+1;
        u_b_out(k,:)=u_b;
        u_f_out(k,:)=u_f;
    end
end

%Plotting
figure;
subplot(2,2,1);       
plot(xx,u_b_out(1,:));
title('Implicit, timestep 2');
axis([0 L 0 2]);
subplot(2,2,2);
plot(xx,u_b_out(2,:));
title('Implicit, timestep 10');
axis([0 L 0 2]);
subplot(2,2,3);
plot(xx,u_b_out(3,:));
title('Implicit, timestep 25');
axis([0 L 0 2]);
subplot(2,2,4);
plot(xx,u_b_out(4,:));
title('Implicit, timestep 100');
axis([0 L 0 2]);
figure;
subplot(2,2,1);       
plot(xx,u_f_out(1,:));
title('Explicit, timestep 2');
axis([0 L 0 2]);
subplot(2,2,2);
plot(xx,u_f_out(2,:));
title('Explicit, timestep 10');
axis([0 L 0 2]);
subplot(2,2,3);
plot(xx,u_f_out(3,:));
title('Explicit, timestep 25');
axis([0 L 0 2]);
subplot(2,2,4);
plot(xx,u_f_out(4,:));
title('Explicit, timestep 100');
axis([0 L 0 2]);

