function [objFun] = th(x,ngrid)
%% PARAMETERS
% INITIAL VALUES
 zst    = x(1); F_0    = x(2); %c_ub   = x(3); 
 q_lb   = x(3);
 % CALL PARAMETERS
param   = parameters()  ; % Call parameters function
% FIXED
r       = param.r       ;   % Monthly discount rate
p       = param.p       ;   % Normalized price
K_s     = param.K_s     ;   % Mg. seller monthly flow cost (Levitt&Venkatesh 2000)
% ESTIMATED 
% Frequency rates
alph    = param.alph    ;   % Rate a buyer meets a seller
gamm    = param.gamm    ;   % Rate of consumption of matched buyer
delt    = param.delt    ;   % Matches' destruction rate
% Consumers' preferences distribution and entry decision
mu_z    = param.mu_z    ;   % Mean of mg. ut distribution
sigz    = param.sigz    ;   % St.dv. of mg. ut distribution
z_ub    = param.z_ub    ;   % Highest drug taste on the support
K_b     = param.K_b     ;   % Mg. buyer entry condition
lamb    = param.lamb    ;   % Fraction of pop. with z = 0 taste for crack
% Consumers' reservation quality distribution
mu_R    = param.mu_R    ;   % Mean of reservation quality
sigR    = param.sigR    ;   % St.dv of reservation quality
R_lb    = param.R_lb    ;   % Lowest reservation quality (from p/z_ub)
R_ub    = p/zst         ;
% Producers' costs distribution and entry
cst     = param.cst     ;   % Mg. seller marginal cost 
xi      = param.xi      ;   % Shape of mg. cost distribution
c_lb    = param.c_lb    ;   % c lower bar
% Selection and measurement error 
mu_et   = param.mu_et   ;   % Selection into ADAM mean  
siget   = param.siget   ;   % Selection into ADAM st.dv.
sigep   = param.sigep   ;   % Quality measurement error st.dv
signu   = param.signu   ;   % Purchases measurement error s.tdv

%% DECLARE DISTRIBUTION FUNCTIONS
% Potential sellers marginal cost
cdfD    = @(x) ((x/cst).^xi).*(x>=c_lb).*(x<=cst)                      ;
pdfD    = @(x) (xi/cst).*((x/cst).^(xi-1)).*(x>=c_lb).*(x<=cst)        ;
% Effective consumers' reservation qualities

cdfM    = @(x) logncdf(x,mu_z,sigz);
%cdfH    = @(x)  (logncdf(x,log(p)-mu_z,sigz)/logncdf(R_ub,log(p)-mu_z,sigz)).*(x<=R_ub)+1*(x>R_ub);            
cdfH    = @(x) (1-cdfM(x))./(1-cdfM(zst)).*(x<=R_ub)+1*(x>R_ub);
%cdfH    = @(x)  (logncdf(x,log(p)-mu_z,sigz)/logncdf(R_ub,log(p)-mu_z,sigz)).*(x<=R_ub)+1*(x>R_ub);            
pdfH    = @(x)  (lognpdf(x,log(p)-mu_z,sigz)/logncdf(R_ub,log(p)-mu_z,sigz)).*(x<=R_ub);
%% DIFFERENTIAL EQUATION - SELLER'S OPTIMAL QUALITY DECISION
% Solve differential equation to find optimal quality q*(c) offered by each seller
cgrid   = ([1:1:ngrid]/ngrid).^(1/xi).*cst ;
cgrid  = (cgrid.')                   ;
%cgrid1  = [.001:(min(cgrid)-.001)/100:min(cgrid),cgrid(2:ngrid)];
%cgrid   = flip(cgrid);
q_0     = q_lb                              ; % Initial condition (for the highest marginal cost seller that enters)
opts    = odeset('RelTol',1e-2,'AbsTol',1e-4);
[c_d,q_d] = ode45(@qdiff,cgrid,q_0,opts,cdfD,pdfD,cdfH,pdfH,F_0,gamm,alph,delt,p);

%% CHECK MONOTONICITY, OTHERWISE CDF CAN'T BE EVALUATED  
qfunc = sort(q_d)   ;
cfunc = c_d         ;
if any(diff(qfunc)<=0) || any(diff(cfunc)>=0) || any(isnan(qfunc)) || any(isinf(qfunc))
    objFun = 10^10;
    return;
end
%%
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

% C. Define f functions (integrands)
fun1    = @(x) x.*pdf_f(x);
fun2    = @(y) phi.*(1 - cdf_F(y))./(r + delta + alpha*(1-cdf_F(y))) ;

%% EQUILIBRIUM CONDITIONS FOR MARGINAL BUYER AND MARGINAL SELLER
% Functions to integrate
f1    = @(x) x.*pdff(x)                                        ;
f2    = @(z) gamm.*(1-cdfF(z))./(r+delt+alph.*(1-cdfF(z)))     ;

% Marginal buyer equation
f(1)  = 100*(alph.*(zst.*(integral(f1,q_lb,q_ub) + integral(f2,p/zst,q_ub) - p)) - K_b)  ;
% Zero profits condition for the marginal seller offering positive quality 
f(2)  = 100*((p-cst.*q_lb).*(1 + (gamm.*delt.*cdfH(q_lb))./(delt + alph.*(1-F_0)))./p - 1);
% FOC for the marginal seller offering positive quality
f(3)  = 100*((p-cst.*q_lb).*((gamm.*delt.*pdfH(q_lb))./(delt + alph.*(1-F_0)).^2) ...
        -cst.*(1 + (gamm.*delt.*cdfH(q_lb))./(delt + alph.*(1-F_0))));
% CDF integrate out to 1
f(4)  = 100*(integral(pdff,q_lb,q_ub)./(1-F_0) - 1)                  ; 

pena    = 10^7*((R_ub<q_lb) + (q_ub<q_lb) + (q_ub<R_ub));

if ~isreal(f)   
    objFun = 10e10;  % Penalty: 
else
    objFun = sum(f.^2) + pena;  % Penalty: 
end

end