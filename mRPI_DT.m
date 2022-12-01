function varargout = mRPI_DT(prob,Phi,K,varargin)
% Solve mRPI subject to discrete-time system
% check the required toolbox
ok_mosek = checkDependency('mosek');
ok_yalmip = checkDependency('yalmip');
if ~ok_yalmip && ~ok_mosek
    error('You need either MOSEK or yalmip installed to use this function.');
end
ok_spotless = checkDependency('spotless');
if ~ok_spotless, error('spotless install is needed'); end


%% ==========begin =============
dimx = 1:prob.dim(1);
dimu = prob.dim(1)+1:prob.dim(1)+prob.dim(2);
dimw = prob.dim(1)+prob.dim(2)+1:sum(prob.dim);
x = prob.var(dimx); 
u = prob.var(dimu);
w = prob.var(dimw);   
%% handle options
if (nargin>3), options=varargin{1}; else, options=struct(); end
if (~isfield(options,'method')), options.method={'SOS'}; end

switch options.methods
    case 'SOS'  
        if (~isfield(options,'max_iterations')), options.max_iterations=100; end
        if (~isfield(options,'converged_tol')), options.converged_tol=.01; end
        if (~isfield(options,'degL')), options.degL = 4; end
        if (~isfield(options,'backoff')), options.backoff = 0; end 
        
        gw = [w-prob.blc(dimw);prob.buc(dimw)-w]; 
        gx = [x-prob.blc(dimx);prob.buc(dimx)-x;]; 
        if prob.feedback
            gu = [-K*x-prob.blc(dimu);prob.buc(dimu)+K*x;]; 
            gx = [gx;gu]; 
        end
        % add feedback
        f = subs(prob.sys,u,-K*x);
        
        [Phi,iterVOL] = RPIsolveSOS(x,w,f,gw,gx,Phi,options); 
        
        varargout{1} = Phi; 
        varargout{2} = iterVOL;  
    case 'SOScontrol'
        if (~isfield(options,'max_iterations')), options.max_iterations=100; end
        if (~isfield(options,'converged_tol')), options.converged_tol=.0001; end
        if (~isfield(options,'degL')), options.degL = 4; end
        if (~isfield(options,'degU')), options.degU = 1; end
        if (~isfield(options,'degLU')), options.degLU = 2; end
        
        gw = [w-prob.blc(dimw);prob.buc(dimw)-w];
        gx = [x-prob.blc(dimx);prob.buc(dimx)-x;];
        gu = [u-prob.blc(dimu);prob.buc(dimu)-u;]; 
        f = prob.sys;
        uf = clean(-K*x);
        [Phi,u]  = mRPI_DT_optcontrol(x,w,u,f,gw,gx,gu,uf,Phi,options); 
        
        varargout{1} = Phi; 
        varargout{2} = u; 
    case 'Minkowski'
        f = prob.sys;
        A = double(diff(f,x));
        B = double(diff(f,u));  
        wuc = prob.buc(dimw); wlc = prob.blc(dimw);
        Wbound = [wuc;wlc];
        n = length(Wbound);
        W_vertex = Wbound;
        for i=1:n/2-1
            W_vertex = [W_vertex,Wbound([i+1:n,1:i])];
        end
%         W_vertex = [0.1, 0.1; 0.1, -0.1; -0.1, -0.1; -0.1, 0.1];
        y = RPIsolveMinkowski(A,B,K,W_vertex); 
        varargout{1} = y;
end
  
end

%% ================ mRPI via SOS ================
function [Phi,iterVOL] = RPIsolveSOS(x,w,f,gw,gx,Phi,options)

Lmonom = monomials([x;w],0:options.degL);
vol = det(Phi);
iterVOL = -log(sqrt(vol)); 
iterV = vol;  
hhh=figure;

for i=1:options.max_iterations
    last_vol = vol;
    [L,rho] = FindL(x,w,f,Phi,gw,gx,Lmonom,options);
    Phi = FindV(x,w,f,L,gw,gx,rho,Phi/rho,Lmonom,options); %*options.backoff 
    vol = det(Phi);
    iterV = [iterV,vol]; 
    % check for convergence
    disp([num2str(vol - last_vol),'<',num2str(options.converged_tol*last_vol )])
    if (vol - last_vol) < options.converged_tol*last_vol 
        delete(hhh)
        break;
    end
    preP = Phi;
    iterVOL = [iterVOL,-log(sqrt(vol))];
    plot(iterVOL,'k-*')
    drawnow
end

end

function [L,rho] = FindL(x,w,f,Phi,gw,gx,Lmonom,options)

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(w);
[prog,rho] = prog.newPos(1);
[prog,L] = prog.newFreePoly(Lmonom);

[prog,Lw] = prog.newFreePoly(Lmonom,length(gw));
[prog,Lx] = prog.newFreePoly(Lmonom,length(gx));

p1 = x'*Phi*x;
p2 = f'*Phi*f;

% setup SOS constraints
prog = prog.withSOS(L);
prog = prog.withSOS(Lw);
prog = prog.withSOS(Lx);
prog = prog.withSOS(rho-p2-L*(1-p1)-Lw'*gw-Lx'*gx);%-gama*p1 
% solve SOS
solver = @spot_mosek;%spot_gurobi;%
pars = spot_sdp_default_options();
pars.verbose = 0;
costfun = rho;
sol = prog.minimize(rho,solver,pars);  %　no optimaize obj

% back off
if options.backoff~=0
    costfun_opt = double(sol.eval(costfun));
    prog = prog.withPos(costfun_opt + (options.backoff)*abs(costfun_opt) - costfun);
    sol = prog.minimize(0,@spot_mosek,pars);
end
L = sol.eval(L); 
rho = double(sol.eval(rho)); 


end


function Phi = FindV(x,w,xp,L,gw,gx,rho,Phi,Lmonom,options)

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(w);

[prog,P] = prog.newPSD(length(x));

p1 = x'*P*x;
p2 = xp'*P*xp;

[prog,Lw] = prog.newFreePoly(Lmonom,length(gw));
[prog,Lx] = prog.newFreePoly(Lmonom,length(gx)); 

% setup SOS constraints
prog = prog.withSOS(Lw);
prog = prog.withSOS(Lx);
prog = prog.withSOS(1-p2-L*(1-p1)-Lw'*gw-Lx'*gx);  
% solve SOS
solver = @spot_mosek;
pars = spot_sdp_default_options();
pars.verbose = 0;
costfun = -trace(Phi^-1*P);
sol = prog.minimize(costfun,solver,pars); %-det(P) replaced by -trace(Phi^-1*P)

% back off
if options.backoff~=0
    costfun_opt = double(sol.eval(costfun));
    prog = prog.withPos(costfun_opt + (options.backoff)*abs(costfun_opt) - costfun);
    sol = prog.minimize(0,@spot_mosek,pars);
end
 
if ~sol.isDualFeasible 
    keyboard
end
if ~sol.isPrimalFeasible
    keyboard
end 
Phi = double(sol.eval(P)); 
end
%% ================ mRPI via SOS control further min ================
function [Phi,uf]  = mRPI_DT_optcontrol(x,w,u,f,gw,gx,gu,uf,Phi,options)
  % ==========begin ===================================
Lmonom = monomials([x;w],0:options.degL);
Lmom = monomials(x,0:options.degL);
um = monomials(x,1:options.degU); 
vol = det(Phi);
iterVOL = -log(sqrt(vol)); 
iterV = vol;  
hhh=figure;

for i=1:options.max_iterations
    last_vol = vol;
    [L,Lu,rho] = FindLU(x,u,w,f,Phi,gw,gx,gu,Lmonom,Lmom,um,uf,options);
    [uf,rho] = FindU(x,u,w,f,L,Lu,Phi,gw,gx,gu,Lmonom,um,uf,options);
    Phi = FindVU(x,u,w,f,L,Lu,gw,gx,Phi/rho,uf,gu,Lmonom,options);
    vol = det(Phi);
    iterV = [iterV,vol]; 
    % check for convergence
    disp([num2str(vol - last_vol),'<',num2str(options.converged_tol*last_vol )])
    if (vol - last_vol) < options.converged_tol*last_vol 
        delete(hhh)
        break;
    end
    preP = Phi;
    iterVOL = [iterVOL,-log(sqrt(vol))];
    plot(iterVOL,'k-*')
    drawnow
    Phi
    uf
end 
end




function [L,Lu,rho] = FindLU(x,u,w,f,Phi,gw,gx,gu,Lmonom,Lmom,um,u0,options)
 
prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(w);

[prog,rho] = prog.newPos(1);
% Create u
uf = [];
for j = 1:size(u,1)
    [prog,un] = prog.newFreePoly(um);
    uf = [uf;un];
end

[prog,L] = prog.newFreePoly(Lmonom);  
[prog,Lw] = prog.newFreePoly(Lmonom,length(gw));
[prog,Lx] = prog.newFreePoly(Lmonom,length(gx)); 
[prog,Lu] = prog.newFreePoly(Lmom,length(gu));

B = options.B;
p1 = x'*Phi*x;
p2 = f'*Phi*f + 2*f'*Phi*B*uf + u0'*B'*Phi*B*uf;  

% setup SOS constraints
prog = prog.withSOS(L); 
prog = prog.withSOS(Lw); 
prog = prog.withSOS(Lx);
prog = prog.withSOS(Lu);
prog = prog.withSOS(1-p2-L*(1-p1)-Lw'*gw-Lx'*gx); 
% prog = prog.withPos(rho-0.9); % max 80% desect
Gu = subs(gu,u,uf);  
prog = prog.withSOS(Gu-Lu*(1-p1)); 
% solve SOS
solver = @spot_mosek;%spot_gurobi;%
pars = spot_sdp_default_options();
pars.verbose = 0;
costfun = 0;
sol = prog.minimize(costfun,solver,pars);  

% back off
% if options.backoff~=0
%     costfun_opt = double(sol.eval(costfun));
%     prog = prog.withPos(costfun_opt + (options.backoff)*abs(costfun_opt) - costfun);
%     sol = prog.minimize(0,@spot_mosek,pars);
% end

L = sol.eval(L);  
Lu = sol.eval(Lu);
rho = double(sol.eval(rho));
% uf = sol.eval(uf);
 
end

function [uf,rho] = FindU(x,u,w,f,L,Lu,Phi,gw,gx,gu,Lmonom,um,u0,options) 
prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(w);

[prog,rho] = prog.newPos(1);
% Create u
uf = [];
for j = 1:size(u,1)
    [prog,un] = prog.newFreePoly(um);
    uf = [uf;un];
end
 
[prog,Lw] = prog.newFreePoly(Lmonom,length(gw));
[prog,Lx] = prog.newFreePoly(Lmonom,length(gx));
  
B = options.B;
p1 = x'*Phi*x;
p2 = f'*Phi*f + 2*f'*Phi*B'*uf + u0'*B'*Phi*B*uf;    

% setup SOS constraints 
prog = prog.withSOS(Lw);
prog = prog.withSOS(Lx);  
prog = prog.withSOS(rho-p2-L*(rho-p1)-Lw'*gw-Lx'*gx); %
prog = prog.withPos(rho-0.8); % max 80% desect
Gu = subs(gu,u,uf);  
prog = prog.withSOS(Gu-Lu*(rho-p1)); 
% solve SOS
solver = @spot_mosek;%spot_gurobi;%
pars = spot_sdp_default_options();
pars.verbose = 0;
costfun = rho;
sol = prog.minimize(costfun,solver,pars); 

% back off
if options.backoff~=0
    costfun_opt = double(sol.eval(costfun));
    prog = prog.withPos(costfun_opt + (options.backoff)*abs(costfun_opt) - costfun);
    sol = prog.minimize(0,@spot_mosek,pars);
end
 
rho = double(sol.eval(rho));
uf = clean(sol.eval(uf));
 
end

function Phi = FindVU(x,u,w,f,L,Lu,gw,gx,Phi,uf,gu,Lmonom,options)
prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(w);

[prog,P] = prog.newPSD(length(x)); 
[prog,Lw] = prog.newFreePoly(Lmonom,length(gw));
[prog,Lx] = prog.newFreePoly(Lmonom,length(gx));

B = options.B;
p1 = x'*P*x;
p2 = f'*P*f + 2*f'*P*B'*uf + uf'*B'*P*B*uf;    

% setup SOS constraints
prog = prog.withSOS(Lw);
prog = prog.withSOS(Lx); 
prog = prog.withSOS(1-p2-L*(1-p1)-Lw'*gw-Lx'*gx );%
Gu = subs(gu,u,uf);  
prog = prog.withSOS(Gu-Lu*(1-p1));
% solve SOS
solver = @spot_mosek;
pars = spot_sdp_default_options();
pars.verbose = 0;  
costfun = trace(Phi^-1*P);
% [prog,tauk] = maxdet(prog,P); 
% costfun = -tauk; 
sol = prog.minimize(costfun,solver,pars);

% back off
if options.backoff~=0
    costfun_opt = double(sol.eval(costfun));
    prog = prog.withPos(costfun_opt + (options.backoff)*abs(costfun_opt) - costfun);
    sol = prog.minimize(0,@spot_mosek,pars);
end

if ~sol.isDualFeasible 
    keyboard
end
if ~sol.isPrimalFeasible
    keyboard
end 
Phi = double(sol.eval(P)); 

end
%% ================ mRPI via Minkowski sum ================
function Z = RPIsolveMinkowski(A,B,K,W_vertex)
 
Ak = A - B*K; 
W = Gao_Polyhedron(W_vertex);
Z = W;
for i = 1:100
    Z1 = Ak^i*W;
    Z = Z+Z1;  
end  
end








