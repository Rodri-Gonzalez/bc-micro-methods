function [objFun] = th(x,alpha,scen)
%% PARAMETERS
% CALL PARAMETERS FUNCTION
param   = parameters()  ; % Call parameters function
% CHANGE NAMES
 zst = x(1); F_0 = x(2); cst = x(3); q_lb = x(4);
% CONDITIONAL ON BASELINE OR COUNTERFACTUAL
if scen == 1 
        alph = alpha; c_ub = x(5); 
elseif scen == 2
        alph = alpha; c_ub = param.c_ub;
end
% FIXED
r       = param.r       ;   % Monthly discount rate
p       = param.p       ;   % Normalized price
K_s     = param.K_s     ;   % Mg. seller monthly flow cost (Levitt&Venkatesh 2000)
S_bar   = param.S_bar   ;   % Potential number of sellers
% ESTIMATED 
% Frequency rates
gamm    = param.gamm    ;   % Rate of consumption of matched buyer
delt    = param.delt    ;   % Matches' destruction rate
% Consumers' preferences distribution and entry decision
B_bar   = param.B_bar   ;
mu_z    = param.mu_z    ;   % Mean of mg. ut distribution
sigz    = param.sigz    ;   % St.dv. of mg. ut distribution
K_b     = param.K_b     ;   % Mg. buyer entry cost
% Consumers' reservation quality distribution
R_ub    = p/zst         ;
% Producers' costs distribution and entry
%cst    = param.cst    ;   % Mg. seller marginal cost 
xi      = param.xi      ;   % Shape of mg. cost distribution
c_lb    = param.c_lb    ;   % c lower bar
%% DECLARE DISTRIBUTION FUNCTIONS
% Effective consumers' reservation qualities
cdfDb   = @(x) ((x/c_ub).^xi).*(x>=c_lb).*(x<=c_ub) + 1.*(x>c_ub)	    ;
cdfD    = @(x) cdfDb(x)./cdfDb(cst).*(x>=c_lb).*(x<=cst) + 1.*(x>cst)   ;
pdfD    = @(x) ((xi/c_ub).*((x/c_ub).^(xi-1)))./(cdfDb(cst)).*(x>=c_lb).*(x<=c_ub)  ;
cdfM    = @(x) logncdf(x,mu_z,sigz);
cdfH    = @(x) (1-cdfM(x))./(1-cdfM(zst)).*(x<=R_ub)+1*(x>R_ub);
pdfH    = @(x)  (lognpdf(x,log(p)-mu_z,sigz)/logncdf(R_ub,log(p)-mu_z,sigz)).*(x<=R_ub);


%% DIFFERENTIAL EQUATION - SELLER'S OPTIMAL QUALITY DECISION
m      = 10000;
cgrid  = ([1:1:m]/m).^(1/xi).*c_ub ;
cgrid  = (cgrid.');
cgrid  = flip(cgrid);
q_0    = q_lb  ; % Initial condition (for the highest marginal cost seller that enters)
opts   = odeset('RelTol',1e-3,'AbsTol',1e-4);
[c_d,q_d] = ode45(@qdiff,cgrid,q_0,opts,cdfD,pdfD,cdfH,pdfH,F_0,gamm,alph,delt,p);
%plot(c_d,q_d)

%% CHECK MONOTONICITY, OTHERWISE CDF CAN'T BE EVALUATED  
qfunc = sort(q_d)   ;
cfunc = c_d         ;
if any(diff(qfunc)<=0) || any(diff(cfunc)>=0) || any(isnan(qfunc)) || any(isinf(qfunc)) || cst>c_ub
    objFun = 10^10;
    return;
end
q       = interp1(cfunc,qfunc,cgrid)        ;
c       = @(q) interp1(qfunc,cfunc,q)       ;
q_ub    = max(q)                            ;   

%% STREET QUALITY DISTRIBUTION (FOR FIRST TIME BUYERS)
% Solve differential equation to find optimal quality q*(c) offered by each seller
cdfF0  = @(x) 1 - cdfD(c(x))*(1-F_0)       ;
cdfF   = @(x) max(F_0,cdfF0(x))            ;  

pdff0  = @(q) ((1-F_0).*pdfD(c(q))./abs(- (2.*gamm.*delt.*(p./c(q) - q).*cdfH(q).*alph.*(1-F_0).*pdfD(c(q)))./...
              ((delt + alph.*(1-F_0).*cdfD(c(q))).^3 + gamm.*delt.*cdfH(q).*(delt + alph.*(1-F_0).*cdfD(c(q))) - ...
               gamm.*delt.*(p./c(q) - q).*pdfH(q).*(delt + alph.*(1-F_0).*cdfD(c(q)))))).*(q>=q_lb).*(q<=q_ub);
           
pdff   = @(q) max(0,pdff0(q))              ;

% B. Define distributions
pdf_f1  = @(q) ((1-F_0).*pdf_d(c(q))./abs(- (2.*phi.*delta.*(p./c(q) - q).*cdf_H(q).*alpha.*(1-F_0).*pdf_d(c(q)))./...
              ((delta + alpha.*(1-F_0).*cdf_D(c(q))).^3 + phi.*delta.*cdf_H(q).*(delta + alpha.*(1-F_0)*cdf_D(c(q))) - ...
               phi.*delta.*(p./c(q) - q).*pdf_h(q).*(delta + alpha.*(1-F_0).*cdf_D(c(q)))))).*(q>=q_lbar) .*(q<=q_ubar);
pdf_f   = @(q) max(0,pdf_f1(q));

% C. Define f functions (integrands) and number of sellers
fun1    = @(x) x.*pdf_f(x);
fun2    = @(y) phi.*(1 - cdf_F(y))./(r + delta + alpha*(1-cdf_F(y))) ;
S       = B_bar*(1-cdfM(zst))*alph/(K_s/p);

%% EQUILIBRIUM CONDITIONS FOR MARGINAL BUYER AND MARGINAL SELLER
% Functions to integrate
f1    = @(x) x.*pdff(x)                                        ;
f2    = @(z) gamm.*(1-cdfF(z))./(r+delt+alph.*(1-cdfF(z)))     ;
% Marginal buyer equation
f(1)  = zst.*(integral(f1,q_lb,q_ub) + integral(f2,p/zst,q_ub))/(p + K_b/alph) - 1  ;  
% Zero profits condition for the marginal seller offering positive quality 
f(2)  = (p-cst.*q_lb).*(1 + (gamm.*delt.*cdfH(q_lb))./(delt + alph.*(1-F_0)))./p - 1;
% FOC for the marginal seller offering positive quality
f(3)  = ((p-cst.*q_lb).*((gamm.*delt.*pdfH(q_lb))./(delt + alph.*(1-F_0)).^2))./ ...
        ((cst.*(1 + (gamm.*delt.*cdfH(q_lb))./(delt + alph.*(1-F_0)).^2)))  - 1;
% CDF integrate out to 1
f(4)  = integral(pdff,q_lb,q_ub)./(1-F_0) - 1                  ; 
% Given an hipotetic number of potential sellers, identify c_ub
f(5)  = S_bar - S*(1-F_0)/cdfDb(cst);
pena  = 10^7*((R_ub<q_lb) + (q_ub<q_lb) + (q_ub<R_ub));

if ~isreal(f)   
    objFun = 10e10;  % Penalty: 
else
    objFun = sum(f.^2) + pena;  % Penalty: 
end

end