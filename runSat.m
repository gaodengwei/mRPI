% test
clear all
close all
clc
dbstop if error
% build system
n = 3; m = 3; dt = 0.2;
x = msspoly('x',n);
u = msspoly('u',m);
w = msspoly('w',n);
x0 = zeros(n,1);
u0 = zeros(m,1);
% Generate a funciton f0 which uses simple arithmetic.

I = diag([5 3 2]);  
fcn = @(omega,T)-inv(I)*cross(omega,I*omega)+ inv(I)*T;
% xn = c2d(@(a)fcn(a,u),x,dt);
xn = c2d(@(a)fcn(a,u*0),x,dt);
xn = xn+fcn(0*x,u)*dt;
AB = double(subs(diff(xn,[x;u]),[x;u],[x0;u0]));
A = AB(1:n,1:n);
B = AB(1:n,n+(1:m));
Q = 10*eye(n);R = eye(m);
[K,S] = dlqr(A,B,Q,R);
uT = -K*x;
xnfd = clean(subs(xn,u,uT));
AB = double(subs(diff(xnfd,[x;u]),[x;u],[x0;u0]));
An = AB(1:n,1:n);
Bn = AB(1:n,n+(1:m));

Phi = 0.0005*S;
V0 = x'*Phi*x;
figure
hold on  
grid on
options.num_samples = 1000;
% y = getLevelSet(x,V0,options);
% plotEllipsoid(y);
% xlabel('x_1')
% ylabel('x_2') 
% zlabel('x_3') 
% axis auto
%% ==============================Minkowski=============================
prob.buc = [20;20;20;10;10;10;0.1;0.2;0.3];
prob.blc = -prob.buc;
prob.dim = [3;3;3];
prob.var = [x;u;w];
prob.sys = An*x+w;   

options.methods = 'Minkowski'; %'Minkowski', 'SOS', 'SOScontrol'
Poly  = mRPI_DT(prob,Phi,0,options);    
[~,volume1] = convexHull(delaunayTriangulation(Poly.V));

% h_2 = plot3(Poly.V(:,1),Poly.V(:,2),Poly.V(:,3),'k:','LineWidth',2);
h_Minkowski = plotEllipsoid(Poly.V',[0,0,0]);
set(h_Minkowski,'FaceAlpha', 0.3)%'EdgeColor','none',
axis([-10,10,-10,10,-10,10])
%% ==============================SOS==========================
% construct the opt problem  for computing mRPI  
clear prob
prob.buc = [20;20;20;10;10;10;0.1;0.2;0.3];
prob.blc =  -prob.buc;
prob.dim = [3;3;3];
prob.var = [x;u;w];
prob.sys = xnfd+w;%
prob.feedback = 0; 
prob.uopt = -K*x; 
options.methods = 'SOS'; 
options.backoff = 0;
options.converged_tol = 0.03;  
% Phi1  = mRPI_DT(prob,Phi,K,options);  
Phi1 = [0.4097    0.0000   -0.0000
    0.0000    0.2769    0.0000
   -0.0000    0.0000    0.1830];
V1 = x'*Phi1*x;
y = getLevelSet(x,V1,options);
dt = delaunayTriangulation(y');
[~,volume2] = convexHull(delaunayTriangulation(y'));

hSOS = plotEllipsoid(y,[1 0 0]);
set(hSOS,'EdgeColor','none','FaceAlpha', 0.3)
xlabel('\omega_1')
ylabel('\omega_2') 
zlabel('\omega_3') 
axis equal
 
keyboard
%% ==============================SOScontrol==========================
options.degU = 2;  
options.methods = 'SOScontrol'; %'Minkowski', 'SOS', 'SOScontrol' 
prob.sys = subs(xn,u,0*u)+w;
options.B = B;
options.backoff = 0.05;
% [Phi3,u]  = mRPI_DT(prob,Phi1,K,options);  
% keyboard
Phi3 = [1.2978    0.0000   -0.0000
    0.0000    0.7867    0.0000
   -0.0000    0.0000    0.7451];
u3 = [ (-11.146)*x(1)+(-0.00036261)*x(2)^2+(-0.00020695)*x(3)^2+(-0.81724)*x(3)*x(2)
        (-0.00033073)*x(1)^2+(-8.6075)*x(2)+(0.00017425)*x(2)*x(1)+(-0.00011454)*x(3)^2+(3.1887)*x(3)*x(1) 
        (-2.431)*x(2)*x(1)+(-8.4564)*x(3)  ];
% [  (-11.146)*x1+(3.2798e-05)*x1^2+(-0.00036261)*x2^2+(3.7615e-05)*x2*x1+(7.617e-05)*x3+(-0.00020695)*x3^2+(-0.81724)*x3*x2  ]
% [  (-0.00033073)*x1^2+(-8.6075)*x2+(-1.3454e-05)*x2^2+(0.00017425)*x2*x1+(1.5973e-05)*x3+(-0.00011454)*x3^2+(3.1887)*x3*x1  ]
% [                      (-5.4946e-05)*x1+(-1.6826e-05)*x2+(-2.431)*x2*x1+(-8.4564)*x3+(7.5176e-05)*x3*x1+(2.7627e-05)*x3*x2  ]
  
  
V3 = x'*Phi3*x;
y = getLevelSet(x,V3,options);
[~,volume3] = convexHull(delaunayTriangulation(y'));

hSOSu = plotEllipsoid(y,[0 0 1]);
set(hSOSu,'EdgeColor','none','FaceAlpha', 1)
xlabel('\omega_1')
ylabel('\omega_2') 
zlabel('\omega_3') 
axis normal
legend('Minkowski','LQR','robust')
set(gcf, 'PaperPosition', [-0.75 0.2 26.5 26]); 
set(gcf, 'PaperSize', [25 25]); 
saveas(gcf, 'WithMargins.pdf'); 
disp(num2str(volume1))
disp(num2str(volume2))
disp(num2str(volume3))


function xn = c2d(f,x,dt)
% rk4 discretization
k1 = f(x)*dt;
k2 = f(x+0.5*k1)*dt;
k3 = f(x+0.5*k2)*dt;
k4 = f(x+0.5*k3)*dt;

xn = clean(x+1/6*(k1+2*k2+2*k3+k4),1e-4);

end







