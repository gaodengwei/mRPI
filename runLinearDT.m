% test
clear all
close all
clc
dbstop if error

A = [1 1;0 1]; 
B = [0.5;1];
Q = eye(2);  
R = 40; 
[K,S] = dlqr(A,B,Q,R); 
% initial moment matrix guess
Phi = 0.1*S;
% build system 
x = msspoly('x',2);
u = msspoly('u',1);
w = msspoly('w',2);
f = A*x+B*u+w;

% plot initial guess of mRPI
figure
hold on  
V0 = x'*Phi*x;
y = getLevelSet(x,V0);
% h_1 = plot(y(1,:),y(2,:),'r--','LineWidth',2);
xlabel('x_1')
ylabel('x_2')
axis auto
% construct the opt problem  for computing mRPI
prob.buc = [5;5;0.3;0.1;0.1];
prob.blc = -prob.buc;
prob.dim = [2;1;2]; % x,u,w
prob.var = [x;u;w];
prob.sys = f; 
prob.uopt = -K*x; 
prob.set0 = [0.1, 0.1; 0.1,-0.1;-0.1,-0.1;-0.1,0.1];


%% ==============================Minkowski=============================
options.methods = 'Minkowski'; %'Minkowski', 'SOS', 'SOScontrol'
Poly  = mRPI_DT(prob,Phi,K,options);    
h_2 = plot(Poly.V(:,1),Poly.V(:,2),'k:','LineWidth',2);
%% ==============================SOS==========================
options.methods = 'SOS'; %'Minkowski', 'SOS', 'SOScontrol'
prob.feedback = 1;
Phi1  = mRPI_DT(prob,0.1*Phi,K,options);  
V1 = x'*Phi1*x;
y = getLevelSet(x,V1);
h_3 = fill(y(1,:),y(2,:),'r');

%% ============================== SOS opt control==========================
options.degU = 2; options.degLU = 4;
options.methods = 'SOScontrol'; %'Minkowski', 'SOS', 'SOScontrol' 
[Phi3,u]  = mRPI_DT(prob,Phi1,K,options);  
V3 = x'*Phi3*x;
y = getLevelSet(x,V3);  
h_4 = fill(y(1,:),y(2,:),'b');
%% =================plot the initial set==========================
y = [-0.1 -0.1;-0.1 0.1;0.1 0.1;0.1 -0.1;-0.1 -0.1]; 
h_0 = plot(y(:,1),y(:,2),'y:','LineWidth',2);
%% =================rerange the coverage==========================
set(gca,'children',[h_0,h_2,h_4,h_3])
legend([h_2,h_0,h_3,h_4],'Minkowski','initial set','LQR' ,'Robust')
 