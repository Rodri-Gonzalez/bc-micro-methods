function [q,c,flag,func] = optimq(c,initial)
zst     = initial(1);
F_0     = initial(2);
cst     = initial(3);
q_lb    = initial(4);

param = parameters();
c_lb    = param.c_lb;
xi      = param.xi;
p       = param.p;
mu_z    = param.mu_z;
sigz    = param.sigz;
R_lb    = param.R_lb;
R_ub    = p/zst;
gamm    = param.gamm;
alph    = param.alph;
delt    = param.delt;
%% DECLARE DISTRIBUTION FUNCTIONS
% Potential sellers marginal cost
cdfD    = @(x) ((x/cst).^xi).*(x>=c_lb).*(x<=cst)                      ;
pdfD    = @(x) (xi/cst).*((x/cst).^(xi-1)).*(x>=c_lb).*(x<=cst)        ;
cdfM    = @(x) logncdf(x,mu_z,sigz);
% Effective consumers' reservation qualities
%cdfH    = @(x)  (logncdf(x,log(p)-mu_z,sigz)/logncdf(R_ub,log(p)-mu_z,sigz)).*(x<=R_ub)+1*(x>R_ub);            
cdfH    = @(x) (1-cdfM(x))./(1-cdfM(zst)).*(x<=R_ub)+1*(x>R_ub);
pdfH    = @(x)  (lognpdf(x,log(p)-mu_z,sigz)/logncdf(R_ub,log(p)-mu_z,sigz)).*(x<=R_ub);

%% SOLVE DIFFERENTIAL EQUATION 
q_0     = q_lb ; % Initial condition (for the highest marginal cost seller that enters)
opts    = odeset('RelTol',1e-2,'AbsTol',1e-4);
c       = flip(c);

[c_d,q_d] = ode45(@qdiff,c,q_0,opts,cdfD,pdfD,cdfH,pdfH,F_0,gamm,alph,delt,p,R_ub);

qfunc = sort(q_d)   ;
cfunc = c_d         ;

if any(diff(qfunc)<=0) || any(diff(cfunc)>=0) || any(isnan(qfunc)) || any(isinf(qfunc))
    flag = 1;
else
    flag = 0;
end

q       = interp1(cfunc,qfunc,c)            ;
c       = @(q) interp1(qfunc,cfunc,q)       ;
c       = c(q)                              ;

func.cdfH = cdfH;
func.pdfH = pdfH;
func.cdfM = cdfM;
func.cdfD = cdfD;
func.pdfD = pdfD;

end