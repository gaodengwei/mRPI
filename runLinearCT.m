clear
close all
clc
dbstop if error

A = [1 1;0 1];
B = [0 0.5;1 0.5]; 
Q = eye(2);  
R = 50*eye(2); 
[K,S] = lqr(A,B,Q,R); 
% Phi = 0.15*S;  
Phi = [9   4
     4   3];
x = msspoly('x',2);
u = msspoly('u',2);
w = msspoly('w',2);
xdot = A*x+B*u+w;

% plot the initial RPI
figure
hold on 
V = x'*Phi*x;
y = getLevelSet(x,V);
% h_1 = plot(y(1,:),y(2,:),'r--','LineWidth',2);
xlabel('x_1')
ylabel('x_2')
 
% build problem
prob.buc = [5;5;0.5;0.5;0.1;0.1];
prob.blc = -prob.buc;
prob.dim = [2;2;2]; % x,u,w
prob.var = [x;u;w];
prob.sys = xdot; 
prob.uopt = -K*x;
prob.gx0 = [0.1-x;x+0.1];

%% ==============================Minkowski=============================
options.methods = 'Aumann'; %'Aumann', 'SOS'
y = mRPI_CT(prob,Phi,K,options); 
h_2 = plot(y(:,1),y(:,2),'k:','LineWidth',2);
%% ==============================SOS=============================
options.methods = 'SOS';
Phi1 = mRPI_CT(prob,Phi,K,options);  
V1 = x'*Phi1*x;
y = getLevelSet(x,V1); 
h_3 = fill(y(1,:),y(2,:),'r');
% Phi1 = [48.3294   18.2734
%    18.2734   15.1237];  
%% ============================== SOS opt control==========================
options.methods = 'SOScontrol';
Phi2 = mRPI_CT(prob,Phi1,K,options);   
% Phi2 = [50.3400   0.0000
%         0.0000   49.6600];   
f2 = x'*Phi2*x;
y = getLevelSet(x,f2); 
h_4 = fill(y(1,:),y(2,:),'b');
%% =================plot the initial set==========================
y = [-0.1 -0.1;-0.1 0.1;0.1 0.1;0.1 -0.1;-0.1 -0.1]; 
h_0 = plot(y(:,1),y(:,2),'y:','LineWidth',2);
%% =================rerange the coverage==========================
set(gca,'children',[h_2 h_0 h_4 h_3 ])
legend([h_2,h_0,h_3,h_4],'Aumann Integrals','Initial set','LQR','Robust')



