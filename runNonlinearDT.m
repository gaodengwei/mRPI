% test
clear all
close all
clc
dbstop if error
% build system 
x = msspoly('x',2); 
u = msspoly('u',1);
w = msspoly('w',2);   
xn = [x(1)+0.01*x(2);
        -0.016*x(1)^3-0.0098*x(1)+x(2)+u]+w;

An = double(subs(diff(xn,x),[x;w],[0;0;0.001;0.001]));   
Bn = double(subs(diff(xn,u),[x;w],[0;0;0.001;0.001]));   
% initial moment matrix guess
% Phi = [7,3; 3,3];  
Phi = eye(2); 
      
figure
hold on    
% V0 = x'*Phi*x;
% y = getLevelSet(x,V0);
% h_1 = plot(y(1,:),y(2,:),'r--','LineWidth',2);
xlabel('\theta')
ylabel('$\dot{\theta}$','Interpreter','latex')
axis auto
% construct the opt problem  for computing mRPI
addNonlinear = 0.016*(1)^3;
%% ==============================Minkowski=============================
prob.buc = [pi;pi;100;0.001;0.001];
prob.blc = [-pi;-pi;-100;-0.001;-0.001]; 
prob.dim = [2;1;2]; % x,u,w
prob.var = [x;u;w];
prob.sys = An*x+w;   
K = 0;
options.methods = 'Minkowski'; %'Minkowski', 'SOS', 'SOScontrol'
Poly  = mRPI_DT(prob,Phi,K,options);    
h_2 = plot(Poly.V(:,1),Poly.V(:,2),'k:','LineWidth',2);

% --------add nonlinear approx-----
prob.buc = [pi/2;pi/2;100;0.001;0.001+addNonlinear];
prob.blc = [-pi/2;-pi/2;-100;-0.001;-0.001-addNonlinear]; 
prob.dim = [2;1;2]; % x,u,w
prob.var = [x;u;w];
prob.sys = An*x+w;   
K = 0;
options.methods = 'Minkowski'; %'Minkowski', 'SOS', 'SOScontrol'
Poly  = mRPI_DT(prob,Phi,K,options);    
h_4 = plot(Poly.V(:,1),Poly.V(:,2),'b:','LineWidth',2);
 
%% ==============================SOS==========================
clear prob
prob.buc = [pi;pi;100;0.001;0.001];
prob.blc =  [-pi;-pi;-100;-0.001;-0.001];
prob.dim = [2;1;2]; % x,u,w
prob.var = [x;u;w];
prob.sys = xn;% xn; An*x+w
prob.feedback = 0;
K=[0,0];
prob.uopt = -K*x; 
options.methods = 'SOS'; 
options.degL = 7;
Phi1  = mRPI_DT(prob,Phi,K,options);  
V1 = x'*Phi1*x;
y = getLevelSet(x,V1);
h_3 = plot(y(1,:),y(2,:),'r','LineWidth',2);

% plot flow
[x1,x2] = meshgrid(-2:0.2:2,-2:0.2:2);  
umax = x1+0.01*x2+0.001; 
vmax = -0.016*x1.^3-0.0098*x1+0.9999*x2+0.001;
umax = (umax-x1)/0.01;
vmax = (vmax-x2)/0.01; 
sumuv = sqrt(umax.^2+vmax.^2);
u1 = umax./sumuv;
v1 = vmax./sumuv;
hhamx = quiver(x1,x2,u1,v1,0.5);
set(hhamx,'Color',[0.5,0.5,0.5],'AutoScaleFactor',0.5) 
umin = x1+0.01*x2-0.001; 
vmin = -0.016*x1.^3-0.0098*x1+0.9999*x2-0.001;
umin = (umin-x1)/0.01;
vmin = (vmin-x2)/0.01; 
sumuv = sqrt(umin.^2+vmin.^2);
u2 = umin./sumuv;
v2 = vmin./sumuv;
hhmin = quiver(x1,x2,u2,v2,0.5);
set(hhmin,'Color',[0.5,0.5,0.5],'AutoScaleFactor',0.5) 
h = sector(x1,x2,u1,v1,u2,v2,0.1); 
legend([h_2,h_3,h_4],'Minkowski-linear','SOS','Minkowski-nonlinear')

