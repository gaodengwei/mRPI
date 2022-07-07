function varargout = mRPI_CT(prob,Phi,K,varargin)
% Solve mRPI subject to continuous-time system
% check the required toolbox
ok_mosek = checkDependency('mosek');
ok_yalmip = checkDependency('yalmip');
if ~ok_yalmip && ~ok_mosek
    error('You need either MOSEK or yalmip installed to use this function.');
end
ok_spotless = checkDependency('spotless');
if ~ok_spotless, error('spotless install is needed'); end

%% ============ begin ============
if (~isfield(prob,'feedback')), prob.feedback = 1;end

dimx = 1:prob.dim(1);
dimu = prob.dim(1)+1:prob.dim(1)+prob.dim(2);
dimw = prob.dim(1)+prob.dim(2)+1:sum(prob.dim);
x = prob.var(dimx);
u = prob.var(dimu);
w = prob.var(dimw); 
if (~isfield(prob,'gx0'))
    gx0 = [x-prob.blc(dimw);prob.buc(dimw)-x];
else
    gx0 = prob.gx0;
end
%% handle options
if (nargin>3), options=varargin{1}; else, options=struct();end
if (~isfield(options,'method')), options.method={'SOS'};end

switch options.methods
    case 'SOS'
        gw = [w-prob.blc(dimw);prob.buc(dimw)-w];
        gx = [x-prob.blc(dimx);prob.buc(dimx)-x;];
        gu = [-K*x-prob.blc(dimu);prob.buc(dimu)+K*x;];
        if prob.feedback
            gx = [gx;gu]; 
        end
            
        % feedback system
        f = subs(prob.sys,u,-K*x);
        
        if (nargin>3),options=varargin{1}; else,options=struct(); end
        if (~isfield(options,'max_iterations')), options.max_iterations=100; end
        if (~isfield(options,'converged_tol')), options.converged_tol=.001; end
        if (~isfield(options,'degL')), options.degL = 4; end
        if (~isfield(options,'degV')), options.degV = 1; end
        if (~isfield(options,'plotVol')), options.plotVol = 0; end
        
        [Phi,iterVOL] = RPIsolveSOS(x,w,f,gw,gx,gx0,Phi,options);
        varargout{1} = Phi; 
        varargout{2} = iterVOL;  
    case 'SOScontrol'
        if (~isfield(options,'max_iterations')), options.max_iterations=100; end
        if (~isfield(options,'converged_tol')), options.converged_tol=.0001; end
        if (~isfield(options,'degL')), options.degL = 4; end
        if (~isfield(options,'degU')), options.degU = 1; end
        if (~isfield(options,'degLU')), options.degLU = 2; end
        if (~isfield(options,'backoff_percent')) options.backoff_percent = 5; end 
        
        gw = [w-prob.blc(dimw);prob.buc(dimw)-w];
        gx = [x-prob.blc(dimx);prob.buc(dimx)-x;];
        gu = [u-prob.blc(dimu);prob.buc(dimu)-u;]; 
        f = prob.sys;
        uf = -K*x;
        [Phi,u]  = mRPI_CT_optcontrol(x,w,u,f,gw,gx,gu,uf,Phi,gx0,options); 
        
        varargout{1} = Phi; 
        varargout{2} = u; 
    case 'Aumann'
        f = prob.sys;
        A = double(diff(f,x));
        B = double(diff(f,u)); 
        wuc = prob.buc(dimw); wlc = prob.blc(dimw);
        Wbound = [wuc;wlc];
        W_vertex = [Wbound,Wbound([2,3,4,1])];
        y = mRPIsolve_Aumann(A,B,K,W_vertex);
        varargout{1} = y;  
end
end
 
function [Phi,iterVOL] = RPIsolveSOS(x,w,f,gw,gx,gx0,Phi,options)

Lmonom = monomials([x;w],0:options.degL);
Lum = monomials(x,0:options.degL);
vol = det(Phi);

% gengerate V 
V = x'*Phi*x; 

iterVOL = -log(sqrt(vol));
for i=1:options.max_iterations
    last_vol = vol;
    L = FindL(x,w,f,V,gw,gx,Lmonom);
    [V,Phi] = FindV(x,w,f,L,gw,gx,gx0,Phi,Lmonom,Lum); 
    vol = det(Phi);
    % check for convergence
    if ((vol - last_vol) < options.converged_tol*last_vol)
        break;
    end
    i
    iterVOL = [iterVOL,-log(sqrt(vol))]; 
end 
end

function L = FindL(x,w,f,V,gw,gx,Lmonom)

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(w);

[prog,rho] = prog.newPos(1);

[prog,L] = prog.newFreePoly(Lmonom);
[prog,Lw] = prog.newFreePoly(Lmonom,length(gw));
[prog,Lx] = prog.newFreePoly(Lmonom,length(gx));

Vdot = diff(V,x)*f;

% setup SOS constraints
% prog = prog.withSOS(L);  % SOS is a suff condition
prog = prog.withSOS(Lw);
prog = prog.withSOS(Lx);
prog = prog.withSOS(-rho*(x'*x)-Vdot-L*(1-V)-Lw'*gw-Lx'*gx);% 
%
% solve SOS
solver = @spot_mosek;%spot_gurobi;%
pars = spot_sdp_default_options();
pars.verbose = 0;
sol = prog.minimize(rho,solver,pars);  %　no optimaize obj
if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end
if ~sol.isPrimalFeasible
    keyboard
end
L = sol.eval(L); 
end

function [V,Phi] = FindV(x,w,f,L,gw,gx,gx0,Phi,Lmonom,Lum)

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(w);

[prog,P] = prog.newPSD(length(x));
V = x'*P*x;

Vdot = diff(V,x)*f;

[prog,Lw] = prog.newFreePoly(Lmonom,length(gw));
[prog,Lx] = prog.newFreePoly(Lmonom,length(gx));
[prog,Lx0] = prog.newFreePoly(Lum,length(gx0));
% setup SOS constraints
prog = prog.withSOS(Lw);
prog = prog.withSOS(Lx);
prog = prog.withSOS(Lx0);
prog = prog.withSOS(-Vdot-L*(1-V)-Lw'*gw-Lx'*gx);
prog = prog.withSOS(1-V-Lx0'*gx0);
% solve SOS
solver = @spot_mosek;
pars = spot_sdp_default_options();
pars.verbose = 0;
sol = prog.minimize(-trace(Phi^-1*P),solver,pars);
% [prog,t] = maxdet(prog,P);
% sol = prog.minimize(-t,solver,pars);

if ~sol.isDualFeasible 
    keyboard
end
if ~sol.isPrimalFeasible
    keyboard
end
Phi = double(sol.eval(P));
V = x'*Phi*x;
end
%% ================ mRPI via SOS control further min ================
function [Phi,u]  = mRPI_CT_optcontrol(x,w,u,f,gw,gx,gu,uf,Phi,gx0,options)
Lmonom = monomials([x;w],0:options.degL);
Lum = monomials(x,0:options.degL);
um = monomials(x,1:options.degU);

vol = det(Phi);
iterVOL = -log(sqrt(vol));
   
for i=1:options.max_iterations
    last_vol = vol;
    [L,Lu] = FindLU(x,u,w,Phi,f,gw,gx,gu,Lmonom,Lum,um); 
    [uf,rho] = FindU(x,u,w,Phi,f,L,Lu,Lmonom,Lum,um,gw,gx,gu,gx0,options);
    Phi = FindVU(x,u,w,Phi/rho,f,L,Lu,gw,gx,gu,gx0,uf,Lmonom,Lum);   
    vol = det(Phi);
    % check for convergence
    if ((vol - last_vol) < options.converged_tol*last_vol)
        break;
    end
    i
    iterVOL = [iterVOL,-log(sqrt(vol))];
end 
end 

function [L,Lu] = FindLU(x,u,w,Phi,f,gw,gx,gu,Lmonom,Lum,um)

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(w);

[prog,gamma] = prog.newPos(1); 
[prog,L] = prog.newFreePoly(Lmonom); 
[prog,Lw] = prog.newFreePoly(Lmonom,length(gw));
[prog,Lx] = prog.newFreePoly(Lmonom,length(gx));
[prog,Lu] = prog.newFreePoly(Lum,length(gu));
% Create u
uf = [];
for j = 1:size(u,1)
    [prog,un] = prog.newFreePoly(um);
    uf = [uf;un];
end

V = x'*Phi*x;
Vdot = diff(V,x)*subs(f,u,uf);
  
Gu = subs(gu,u,uf);  
% setup SOS constraints 
prog = prog.withSOS(Lw);
prog = prog.withSOS(Lx);
prog = prog.withSOS(Lu);
prog = prog.withSOS(-gamma*(x'*x)-Vdot-L*(1-V)-Lw'*gw-Lx'*gx); 
prog = prog.withSOS(Gu-Lu*(1-V)); 
 
% solve SOS
solver = @spot_mosek;%spot_gurobi;%
pars = spot_sdp_default_options();
pars.verbose = 0; 

sol = prog.minimize(gamma,solver,pars);
if ~sol.isDualFeasible 
    keyboard
end
if ~sol.isPrimalFeasible
    keyboard
end
L = sol.eval(L); 
Lu = sol.eval(Lu);   
end

function [uf,rho] = FindU(x,u,w,Phi,f,L,Lu,Lmonom,Lum,um,gw,gx,gu,gx0,options) 
prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(w);

[prog,rho] = prog.newPos(1);
% [prog,gamma] = prog.newPos(1); 
[prog,Lw] = prog.newFreePoly(Lmonom,length(gw));
[prog,Lx] = prog.newFreePoly(Lmonom,length(gx));
[prog,Lx0] = prog.newFreePoly(Lum,length(gx0));

% Create u
uf = [];
for j = 1:size(u,1)
    [prog,un] = prog.newFreePoly(um);
    uf = [uf;un];
end
  
V = x'*Phi*x;
Vdot = diff(V,x)*subs(f,u,uf);
 
Gu = subs(gu,u,uf);  
% setup SOS constraints
prog = prog.withSOS(Lw);
prog = prog.withSOS(Lx);
prog = prog.withSOS(Lx0);
prog = prog.withSOS(-Vdot-L*(rho-V)-Lw'*gw-Lx'*gx);%-gamma*(x'*x)
prog = prog.withSOS(rho-V-Lx0'*gx0);
prog = prog.withSOS(Gu-Lu*(rho-V)); 
 
% solve SOS
solver = @spot_mosek;
pars = spot_sdp_default_options();
pars.verbose = 0; 
sol = prog.minimize(rho,solver,pars);
if ~sol.isDualFeasible 
    keyboard
end
if ~sol.isPrimalFeasible
    keyboard
end

%% Back off on objective
rhof = double(sol.eval(rho)); 
[prog,bslack] = prog.newPos(1); 
prog = prog.withEqs(bslack - (-(1-options.backoff_percent/100)*rhof + rho));

sol = prog.minimize(0,solver,pars); 
  
uf = sol.eval(uf);  
rho = double(sol.eval(rho)); 
end
 
function  Phi  = FindVU(x,u,w,Phi,f,L,Lu,gw,gx,gu,gx0,uf,Lmonom,Lum)  
prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(w);
 
[prog,P] = prog.newPSD(2); 
V = x'*P*x; 
Vdot = diff(V,x)*subs(f,u,uf);

[prog,Lw] = prog.newFreePoly(Lmonom,length(gw));
[prog,Lx] = prog.newFreePoly(Lmonom,length(gx));
[prog,Lx0] = prog.newFreePoly(Lum,length(gx0));

% setup SOS constraints
prog = prog.withSOS(Lw);
prog = prog.withSOS(Lx);
prog = prog.withSOS(Lx0);
prog = prog.withSOS(-Vdot-L*(1-V)-Lw'*gw-Lx'*gx); %-gamma*(x'*x)
Gu = subs(gu,u,uf);  
prog = prog.withSOS(Gu-Lu*(1-V)); 
prog = prog.withSOS(1-V-Lx0'*gx0);

% solve SOS
solver = @spot_mosek;
pars = spot_sdp_default_options();
pars.verbose = 0;
sol = prog.minimize(-trace(Phi^-1*P),solver,pars); %-det(P) replaced by -trace(Phi^-1*P) 
if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end
 
Phi = double(sol.eval(P)); 
end
%% ================ mRPI via Aumann int ================
function y = mRPIsolve_Aumann(A,B,K,W_vertex) 
Ak = A - B*K;  
syms t
fun = expm(Ak*t)*W_vertex';
Sinf  = W_vertex;
for i=1:20
    Sinf = [Sinf;double(int(fun, t, 0, i))'];
end
chll = convhull(Sinf);
y = Sinf(chll,:); 
end