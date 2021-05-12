function [q,c,flag,func] = optimq(c,initial,alpha,c_ubar)
zst     = initial(1);
F_0     = initial(2);
cst     = initial(3);
q_lb    = initial(4);
alph    = alpha;
c_ub    = c_ubar;

param = parameters();
c_lb    = param.c_lb;
xi      = param.xi;
p       = param.p;
mu_z    = param.mu_z;
sigz    = param.sigz;
R_lb    = param.R_lb;
R_ub    = p/zst;
gamm    = param.gamm;
delt    = param.delt;
%% DECLARE DISTRIBUTION FUNCTIONS

cdfD    = @(x) ((x/cst).^xi).*(x>=c_lb).*(x<=cst)                      ;
pdfD    = @(x) (xi/cst).*((x/cst).^(xi-1)).*(x>=c_lb).*(x<=cst)        ;
cdfM    = @(x) logncdf(x,mu_z,sigz);
          
cdfH    = @(x) (1-cdfM(x))./(1-cdfM(zst)).*(x<=R_ub)+1*(x>R_ub);
pdfH    = @(x) (lognpdf(x,log(p)-mu_z,sigz)/logncdf(R_ub,log(p)-mu_z,sigz)).*(x<=R_ub);

%% SOLVE DIFFERENTIAL EQUATION 
q_0     = q_lb ; % Initial condition (for the highest marginal cost seller that enters)
opts    = odeset('RelTol',1e-3,'AbsTol',1e-4);
c       = flip(c);
[c_d,q_d] = ode45(@qdiff,c,q_0,opts,cdfD,pdfD,cdfH,pdfH,F_0,gamm,alph,delt,p,R_ub);
% Keep unique values
qfunc = sort(q_d)   ;
dq    = diff(qfunc) ;
dq    = [0;dq]      ;
cfunc = c_d         ;
qfunc(dq(:,1)<=0,:)=[];
cfunc(dq(:,1)<=0,:)=[];

if any(diff(qfunc)<=0) || any(diff(cfunc)>=0) || any(isnan(qfunc)) || any(isinf(qfunc))
    flag = 1;
else
    flag = 0;
end

q = interp1(cfunc,qfunc,c)            ;
c = @(q) interp1(qfunc,cfunc,q)       ;
qf= @(c) interp1(cfunc,qfunc,c)       ;

func.cdfH = cdfH;
func.pdfH = pdfH;
func.cdfM = cdfM;
func.cdfD = cdfD;
func.pdfD = pdfD;
func.c    = c;
func.qf   = qf;

end