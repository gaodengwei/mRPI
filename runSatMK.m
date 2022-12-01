clc
clear
close all

fx = @(x1,x2,x3,u1,u2,u3)[(0.04)*u1+x1+(0.0033333)*x2^2*x1+(-0.0033333)*x3^2*x1+(0.04)*x3*x2+(-0.00073333)*x3*x2*x1^2
    (0.066667)*u2+x2+(-0.016667)*x2*x1^2+(-0.2)*x3*x1+(0.001)*x3*x1^3+(0.0002)*x3^3*x1+(-0.0033333)*x3^2*x2+(-0.00073333)*x3*x2^2*x1+(0.00010667)*x3^2*x2*x1^2
    (0.2)*x2*x1+(-0.001)*x2*x1^3+(0.0002)*x2^3*x1+(0.1)*u3+x3+(-0.016667)*x3*x1^2+(0.0033333)*x3*x2^2+(-0.00010667)*x3*x2^2*x1^2+(-0.00073333)*x3^2*x2*x1  ];
u_lqr = @(x)-[2.9686   -0.0000   -0.0000
    -0.0000    2.8465    0.0000
    0.0000   -0.0000    2.7016]*x;
u_sos =  @(x)[ (-11.146)*x(1)+(-0.00036261)*x(2)^2+(-0.00020695)*x(3)^2+(-0.81724)*x(3)*x(2)
        (-0.00033073)*x(1)^2+(-8.6075)*x(2)+(0.00017425)*x(2)*x(1)+(-0.00011454)*x(3)^2+(3.1887)*x(3)*x(1) 
        (-2.431)*x(2)*x(1)+(-8.4564)*x(3)  ];
% u_sos =  @(x)[ (-11.146)*x(1)+(-0.81724)*x(3)*x(2)
%         (-8.6075)*x(2)+(3.1887)*x(3)*x(1) 
%         (-2.431)*x(2)*x(1)+(-8.4564)*x(3)  ];
N = 50;
figure
h1 = subplot(3,1,1);hold on;xlabel('t(s)');ylabel('\omega_1');
h2 = subplot(3,1,2);hold on;xlabel('t(s)');ylabel('\omega_2');
h3 = subplot(3,1,3);hold on;xlabel('t(s)');ylabel('\omega_3');

figure
hh1 = subplot(3,1,1);hold on;
hh2 = subplot(3,1,2);hold on;
hh3 = subplot(3,1,3);hold on;
% figure
% g1 = subplot(3,1,1);hold on;
% g2 = subplot(3,1,2);hold on;
% g3 = subplot(3,1,3);hold on;
figure
gg1 = subplot(3,1,1);hold on;
gg2 = subplot(3,1,2);hold on;
gg3 = subplot(3,1,3);hold on;
 
rm_lqr = 0; rm_sos = 0;
for j = 1:200
    x = 5-10*rand(3,1);
    w = [0.1;0.2;0.3] - 2*[0.1;0.2;0.3].*rand(3,1);
    y = x;
    xn_lqr = x;xn_sos = x;
    un_lqr = []; un_sos = [];
    for i=1:N
        ux = u_lqr(x);ux = min(max(ux,-10),10);
        x = fx(x(1),x(2),x(3),ux(1),ux(2),ux(3))+w;
        uy = u_sos(y);uy = min(max(uy,-10),10);
        y = fx(y(1),y(2),y(3),uy(1),uy(2),uy(3))+w;
        xn_lqr = [xn_lqr,x];un_lqr= [un_lqr,ux];
        xn_sos = [xn_sos,y];un_sos= [un_sos,uy];
    end
    if ~any(isnan(x(:)))
        plot3(h1,(0:N)*0.2,xn_lqr(1,:),-1*ones(N+1,1),'color',[0.85, 0.33, 0.10])
        plot3(h2,(0:N)*0.2,xn_lqr(2,:),-1*ones(N+1,1),'color',[0.85, 0.33, 0.10])
        plot3(h3,(0:N)*0.2,xn_lqr(3,:),-1*ones(N+1,1),'color',[0.85, 0.33, 0.10])
         
        plot(hh1,(1:N)*0.2,un_lqr(1,:),'r')
        plot(hh2,(1:N)*0.2,un_lqr(2,:),'r')
        plot(hh3,(1:N)*0.2,un_lqr(3,:),'r')
    else
        rm_sos=rm_sos+1;
    end
    if ~any(isnan(y(:)))
        plot3(h1,(0:N)*0.2,xn_sos(1,:),zeros(N+1,1),'color',[0.00, 0.45, 0.74])
        plot3(h2,(0:N)*0.2,xn_sos(2,:),zeros(N+1,1),'color',[0.00, 0.45, 0.74])
        plot3(h3,(0:N)*0.2,xn_sos(3,:),zeros(N+1,1),'color',[0.00, 0.45, 0.74])
        
        plot(gg1,(1:N)*0.2,un_sos(1,:),'b')
        plot(gg2,(1:N)*0.2,un_sos(2,:),'b')
        plot(gg3,(1:N)*0.2,un_sos(3,:),'b')
    else
        rm_lqr=rm_lqr+1;
    end
end
legend('LQR','robust')
keyboard