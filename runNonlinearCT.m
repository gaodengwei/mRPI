% test
clear all
close all
clc
dbstop if error
% build system 
x = msspoly('x',2); 
u = msspoly('u',1);
w = msspoly('w',2);  
xdot = [-x(2)*w(1);
        -x(2)*(1-x(1)^2)+x(1)+w(2)+u];

A = double(subs(diff(xdot,x),[x;w],[0;0;1;0.2]));   
% initial moment matrix guess  
Phi = [1.5   -0.2
   -0.2    1.5]; 
% -------------------- find the bound of disturbance-------------------------
% W = subs(A*x+w-xdot,u,0);
% lim = 0.2;
% [x1,x2,x3,x4] = ndgrid(-0.36:0.2:0.36,-0.36:0.2:0.36,0.9:0.1:1,-0.2:0.2:0.2);
% aaa = msubs(W(1),[x;w],[x1(:),x2(:),x3(:),x4(:)]');
% bbb = msubs(W(2),[x;w],[x1(:),x2(:),x3(:),x4(:)]'); 
% figure;subplot(2,1,1),plot(aaa);subplot(2,1,2),plot(bbb);
% ------------------------------------------------------------


figure
hold on    
V0 = x'*Phi*x;
% y = getLevelSet(x,V0);
% h_1 = plot(y(1,:),y(2,:),'r--','LineWidth',2);
xlabel('x_1')
ylabel('x_2') 
axis auto
% construct the opt problem  for computing mRPI
%% ==============================Aumann=============================
prob.buc = [pi;2*pi;100;1;0.47];
prob.blc = [-pi;-2*pi;-100;-0.87;-0.32]; 
prob.dim = [2;1;2]; % x,u,w
prob.var = [x;u;w];
prob.sys = A*x+w;   
K = 0;
options.methods = 'Aumann'; %'Aumann', 'SOS', 'SOScontrol'
Poly  = mRPI_CT(prob,Phi,K,options);    
h_2 = plot(Poly(:,1),Poly(:,2),'k:','LineWidth',2);
 
% keyboard
%% ==============================SOS==========================
clear prob
prob.buc = [pi;2*pi;100;1.1;0.2];
prob.blc =  [-pi;-2*pi;-100;0.9;-0.2];
prob.dim = [2;1;2]; % x,u,w
prob.var = [x;u;w];
prob.sys = xdot; 
prob.feedback = 0;
K=[0,0];
prob.uopt = -K*x; 
prob.gx0 = [0.01+x(1); 0.01+x(2);0.01-x(1);0.01-x(2)]; 
options.methods = 'SOS'; 
Phi1  = mRPI_CT(prob,Phi,K,options);  
V1 = x'*Phi1*x;
y = getLevelSet(x,V1);
h_3 = plot(y(1,:),y(2,:),'r','LineWidth',2);
 
% plot flow
[x1,x2] = meshgrid(-2:0.2:2,-1.5:0.3:1.5); 
umax = -x2*1.1;
vmax = -x2.*(1-x1.^2)+x1+0.2; 
sumuv = sqrt(umax.^2+vmax.^2);
u1 = umax./sumuv;
v1 = vmax./sumuv;
hhamx = quiver(x1,x2,u1,v1,0.5);
set(hhamx,'Color',[0.5,0.5,0.5],'AutoScaleFactor',0.5)

umin = -x2*0.9;
vmin = -x2.*(1-x1.^2)+x1-0.2; 
sumuv = sqrt(umin.^2+vmin.^2);
u2 = umin./sumuv;
v2 = vmin./sumuv;
hhmin = quiver(x1,x2,u2,v2,0.5);
set(hhmin,'Color',[0.5,0.5,0.5],'AutoScaleFactor',0.5) 
h = sector(x1,x2,u1,v1,u2,v2,0.165); 
legend([h_2,h_3],'Aumann','SOS')

