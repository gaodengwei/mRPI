function g = TransFormSystem(xdot,x,x0)
% from xdot = f(t,x,u) to zdot = g(t,z,u) via z=c(t,x) and x=d(t,z)
 % zdot = dcdt(t,d(t,z)) + dcdx(t,d(t,z)*f(t,d(t,z),u)

 % x0 is the nominal system
c = x - x0;
d = x + x0; 

g = diff(c,x)*xdot;
g = subss(g,x,d);

end

