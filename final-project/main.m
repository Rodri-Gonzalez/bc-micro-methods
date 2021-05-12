%**************************************************************************
%**************************************************************************
% EMPIRICAL METHODS - FINAL PROJECT 
% Simulation of Galenianos & Gavazza (AER 2017) 
% "A Structural Model of the Retail Market for Illicit Drugs" + 
% Matching efficiency counterfacual
% Rodrigo Gonzalez
%**************************************************************************
% The file main.m runs the process. The function th.m estimate the
% thresholds and optimq.m uses the thresholds to estimate the solution to 
% the differential equation that solves the seller decision with the file 
% qdiff.m. Additionally, parameters.m store the parameters of the model 
% extracted from the paper. To run baseline and counterfactual scenario
% just change "scen" on line 25 and "coeff" in line 34
%**************************************************************************
clear all 
cd 'D:\GitDir\bc-micro-methods\final-project'
rand('seed',1116);      % For reproductibility

%% 0) SELECT METHOD AND SCENARIO. IN CASE OF COUNTERFACTUAL, SELECT INCREASE IN OMEGA IN ``COEFF" (line 27)
% Method -- (1 for Fminsearch)
meth = 1;
% Scenario -- (1 for Baseline || 2 for Counterfactual) Change
% counterfactual coefficient in line 34.
scen = 1;
% Call parameters
param   = parameters() ;  % Call parameters
% 
if scen == 1
    % Use meetings rate from the paper
    alph = param.alph;
elseif scen == 2
    % Increase matching efficiency 25%
    coeff = 1.25;
    alph = (param.ome*coeff)^2*param.p/param.K_s;
end
% Call function to estimate thresholds
thresh  = @(x) th(x,alph,scen);

%% 1) ESTIMATE EQUILIBRIUM THRESHOLDS IMPLIED BY THE MODEL [z*,F_0,q_lb,c*,c_ub]
    % Baseline, FMINSEARCH 
if meth == 1 && scen==1
    initial  = [150,0.17,122,0.63,126];
    options  = optimset('PlotFcns','optimplotfval','MaxFunEvals',200000,'MaxIter',300);
    thre_est = fminsearch(thresh,initial,options);
    c_ub    = thre_est(5);
    % Counterfactual, FMINSEARCH
elseif meth == 1 && scen==2   
    options= optimset('PlotFcns','optimplotfval','MaxFunEvals',200000,'MaxIter',300);
    initial = [150,0.17,122,0.63];
    thre_est= fminsearch(thresh,initial,options);
    c_ub    = param.c_ub;
end
% Use estimated parameters
zst = thre_est(1); F_0  = thre_est(2);
cst = thre_est(3); q_lb = thre_est(4); 
theta = [thre_est(1), thre_est(2), thre_est(3), thre_est(4)];

%% 1.5) CALL PARAMETERS STORED IN FUNCTION PARAMETERS()
ngrid   = 5000                           ;  % Number of simulations
% BUYERS
mu_z    = param.mu_z                     ;  % Mean valuation
sigz    = param.sigz                     ;  % Stdv of valuation
B_max   = param.B_max                    ;  % Eq number of consumers
B_bar   = param.B_bar                    ;   
p       = param.p                        ;  % Given price
z       = exp(mu_z + sigz*randn(ngrid,1));  % Random draws of consumers tastes
R       = p./z                           ;  % Reservation qualities                                  
R_lb    = param.R_lb                     ;  % R lower bar (p/z_ub)
% SELLERS
xi      = param.xi                       ;  % c distribution shape parameter
K_s     = param.K_s                      ;  % Sellers entry fixed cost
% RATES
gamm    = param.gamm                     ;
delt    = param.delt                     ; 
% OTHERS
sigep   = param.sigep                    ;
mu_eps  = param.mu_eps                   ;
mu_eta  = param.mu_et                    ;
sigeta  = param.siget                    ;
lamb    = param.lamb                     ;
signu   = param.signu                    ;

%% 2) USE ESTIMATED THRESHOLDS TO SOLVE DIFFERENTIAL EQUATION
cgrid   = ([1:1:ngrid]/ngrid).^(1/xi).*c_ub ;
cgrid   = (cgrid.'); 
[q,c,flag,funct] = optimq(cgrid,theta,alph,c_ub);  
q_ub = max(q);
cdfH = funct.cdfH;
pdfH = funct.pdfH;
cdfM = funct.cdfM;
cdfD = funct.cdfD;
pdfD = funct.pdfD;
c    = funct.c;  
qf   = funct.qf; 
%plot(c(q),q)

%% 3) STREET QUALITY DISTRIBUTION (FOR FIRST TIME BUYERS)
% Plug in differential solution of optimal quality q*(c) offered by each seller
cdf_F = @(q) 1-(1-F_0)*cdfD(c(q));
cdfF  = @(q) max(F_0,cdf_F(q))   ;
% Compute inverse CDF to use in simulation
invF  = @(x) max((x<=F_0).*0, (x>F_0).*(qf(((1-x)./(1-F_0)).^(1/xi).*cst))); 
% Compute pdf of F(q) to calculate moment conditions, use derivative of
% implicit function
pdff0 = @(q) ((1-F_0).*pdfD(c(q))./abs(- (2.*gamm.*delt.*(p./c(q) - q).*cdfH(q).*alph.*(1-F_0).*pdfD(c(q)))./...
              ((delt + alph.*(1-F_0).*cdfD(c(q))).^3 + gamm.*delt.*cdfH(q).*(delt + alph.*(1-F_0).*cdfD(c(q))) - ...
               gamm.*delt.*(p./c(q) - q).*pdfH(q).*(delt + alph.*(1-F_0).*cdfD(c(q)))))).*(q>=q_lb).*(q<=q_ub);        
pdff  = @(q) max(0,pdff0(q));

%% 4) ADD MEASUREMENT ERROR (ONLY ON THE CONTINUOUS PART OF THE QUALITY DISTRIBUTION)            
pdfEps= @(x) lognpdf(x,mu_eps,sigep);
pdfQr = @(x,y) (1./x).*pdfEps(x).*pdff(y./x);
% Integration bounds on epsilon
qlowest = logninv(.00001,mu_eps,sigep); qupper = logninv(.99999,mu_eps,sigep);
% Integration bounds on q
qmin  = q_lb*qlowest; qmax = q_ub*qupper;

%% 5) CORRECT FOR ADAM (ARRESTED USERS) SELECTION TO CALCULATE STATISTICS
num   = @(x) normpdf(x,mu_z,sigz).*(1-normcdf(-x,mu_eta,sigeta));
numb  = @(x) normpdf(x,mu_z,sigz).*normcdf(-x,mu_eta,sigeta);
den   = integral(num,log(zst),mu_z + 5*sigz); 
denb  = integral(numb,log(zst),mu_z+5*sigz);

%% 6) SIMULATE MEETINGS, BEGIN FROM UNMATCHED STATE
% #of consumers searching
cons  = 1000;
% #of periods to simulate
T     = 1000; 
m     = zeros(cons,1);
% Create vectors to store simulation data
z_sim = zeros(cons,1); R_sim = zeros(cons,1);
purch = zeros(cons,1); match = zeros(cons,1);
reg   = zeros(cons,1); q_buy = zeros(cons,1);

for j = 1:cons
    fprintf('%2.4f\tSimulating\n',j)
    
    % Draw random preference for a potential consumer consumer
    %[~,zrnd]  = min(abs(cdf_lnzg-rand(1)));
    %z_sim(j,1)= z(zrnd,1);
     z_sim(j,1)= exp(mu_z + sigz*randn(1));
     R_sim(j,1)= p/z_sim(j,1);
    
     t = 1; % Restart counter for each consumer
     % BEGIN WITH BURN IN PERIOD
        while t < T
        % Consumer ENTERS THE MARKET for the first time UNMATCHED, with Q=0               
        if t>0 && t<(T-30) && match(j,1) == 0
            % UNMATCHED CONSUMERS - BURN IN
            t = t + exprnd(30/alph);
                q_buy(j,1) = invF(rand(1));
                if q_buy(j,1)>=R_sim(j,1)
                    match(j,1) = 1;
                else
                    match(j,1) = 0; % If didn't matched, waits until finding new seller.
                end
        elseif t>0 && t<(T-30) && match(j,1) == 1
            % UNMATCHED CONSUMERS - BURN IN PHASE
            t1 = exprnd(30/alph); t2 = exprnd(30/gamm); t3 = exprnd(30/delt);    
            t  = t + min([t1,t2,t3]); % New event for matched consumers   
            % Meet a new seller
            if min([t1,t2,t3])==t1
                q_off = invF(rand(1));
                    if q_off>=q_buy(j,1)
                    match(j,1) = 1; q_buy(j,1) = q_off; reg(j,1) = 0;
                    else
                    match(j,1) = 1;
                    end
            % Consumes from previous seller         
            elseif min([t1,t2,t3])==t2
                reg(j,1) = 1;    
            elseif min([t1,t2,t3])==t3
                match(j,1) = 0; q_buy(j,1)   = 0; reg(j,1) = 0;   
            end
        
        % BURN IN PHASE ENDS - STARTS RECORDING CONSUMPTION FREQUENCY    
        % FINAL 30 DAYS UNMATCHED
        elseif t>(T-30) && t<T && match(j,1) == 0
            t = t + exprnd(30/alph);
            q_buy(j,1)   = invF(rand(1));
            purch(j,1) = 1 + purch(j,1); % Record purchases
                if q_buy(j,1)>=R_sim(j,1)
                    match(j,1) = 1;
                else
                    match(j,1) = 0; % If didn't matched, waits until finding new seller.
                end   
        % FINAL 30 DAYS MATCHED    
        elseif t>(T-30) && t<T && match(j,1) == 1
            % MATCHED CONSUMERS
            t1 = exprnd(30/alph); t2 = exprnd(30/gamm); t3 = exprnd(30/delt);    
            t  = t + min([t1,t2,t3]); % New event for matched consumers   
            % Meet a new seller
            if min([t1,t2,t3])==t1
                q_off = invF(rand(1));
                purch(j,1) = 1 + purch(j,1); % Record purchases
                    if q_off>=q_buy(j,1)
                    match(j,1) = 1; q_buy(j,1) = q_off; reg(j,1) = 0;
                    else
                    match(j,1) = 1;
                    end
            % Consumes from previous seller         
            elseif min([t1,t2,t3])==t2
                purch(j,1) = 1 + purch(j,1); % Record purchases
                reg(j,1) = 1;    
            elseif min([t1,t2,t3])==t3
                match(j,1) = 0; q_buy(j,1)   = 0; reg(j,1) = 0;   
            end  
        end
        end
end

%% 7) CALCULATE CONSUMPTION STATISTICS
% Arrest rate and people consuming drugs on NSDUH (Health Survey)
arrest = (lamb + (1 - lamb)* cdfM(zst))*(1 - normcdf(0,mu_eta,sigeta)) + (1-lamb)*den;
nsduh  = (1-lamb)*denb/(1-arrest);
purch1 = purch;
% Add purchases measurement error
purch  = round(purch1.*lognrnd(-0.5*signu^2,signu,cons,1));
f30d   = (1-lamb)*den*(1-mean(match==0))/arrest; 
lastreg = sum((reg==1))/sum((purch>0));

meanp_m = sum(purch.*(reg==1).*(purch>0))/sum((reg==1).*(purch>0));
meanp_u = sum(purch.*(reg==0).*(purch>0))/sum((reg==0).*(purch>0));

meanq_m = sum(q_buy.*(reg==1).*(purch>0))/sum((reg==1).*(purch>0));
meanq_u = sum(q_buy.*(reg==0).*(purch>0))/sum((reg==0).*(purch>0));
meanq   = mean(q_buy);

std_m   = std(purch((reg==1)&(purch>0)));
std_u   = std(purch((reg==0)&(purch>0)));

med_m   = median(purch((reg==1)&(purch>0)));
med_u   = median(purch((reg==0)&(purch>0)));

%% 8) CALCULATE STATISTICS OF THE QUALITY DISTRIBUTION
% Fraction of scams/rip-offs
m1 = F_0;
% Mean (Conditional on q>0, includes measurement error)
fm2= @(eps,q) q.*pdfQr(eps,q)/(1-F_0);
m2 = integral2(fm2,qlowest,qupper,qmin,qmax);

% Standard deviation (Conditional on q>0)
fm3= @(eps,q) ((q-m2).^2).*pdfQr(eps,q)/(1-F_0);
m3 = sqrt(integral2(fm3,qlowest,qupper,qmin,qmax));

% Slow median calculation, improve it
% Median (Conditional on q>0, includes measurement error)
%q_med = qmax/12; i=0;
%while i<0.5
%    q_med = q_med + 0.01; 
%    med = integral2(pdfQr,qlowest,qupper,qmin,q_med)/(1-F_0);
%    i = med;
%end
%m4 = q_med;

m4=0;

% Skewness (Conditional on q>0, includes measurement error)
fm5= @(eps,q) (((q-m2)/m3).^3).*pdfQr(eps,q)/(1-F_0);
m5 = integral2(fm5,qlowest,qupper,qmin,qmax);

% Kurtosis (Conditional on q>0, includes measurement error)
fm6= @(eps,q) (((q-m2)/m3).^4).*pdfQr(eps,q)/(1-F_0);
m6 = integral2(fm6,qlowest,qupper,qmin,qmax);

%% 9) CALIBRATE COBB-DOUGLAS MATCHING FUNCTION PARAMETERS AND COMPUTE RELEVANT STATISTICS
% Buyers and sellers
B   = B_max*(1-cdfM(zst));
ome = (alph*K_s/p)^(1/2); 
S   = B*alph^2/ome^2;

% Active-unmatched consumers distribution
fun_un = @(z) B*delt.*pdfH(z)./(delt + alph.*(1 - cdfF(z)));                  
n_r  = @(y) integral(fun_un,0.001,y);        
n_ub = n_r(p/zst);

% Rate of matched consumers
matched_b = (B-n_ub)/B;
trans_0q = alph*B*F_0/(alph*B + B - n_ub)*gamm;

% Quality of matched consumers
pdfG = @(q) alph.*B.*delt.*pdff(q).*cdfH(q)./(((delt + alph.*(1-cdfF(q))).^2).*(B-n_ub));
fun8 = @(q) q.*pdfG(q);
meanq_match = integral(fun8,q_lb,q_ub);
% Mean quality in the market and standard deviation
fm2_s   = @(q) q.*pdff(q); m2_s = integral(fm2_s,qmin,qmax);
m3_s    = std(q_buy(purch1>0));
meanp_t = mean(purch1);
meanq_t = mean(q_buy(purch1>0))*mean(purch1(purch1>0));

%% 10) DISPLAY SIMULATION RESULTS - TABLE 3 OF THE PAPER
format short
disp('==================================================');
disp('Simulated moments                         ');
disp(['Rip-offs (%)                             ', num2str(F_0)]);
disp(['Avg. Quality per $100, (q>0)             ', num2str(m2)]);
disp(['Std.dev. Quality per $100, (q>0)         ', num2str(m3)]);
disp(['Median Quality per $100, (q>0)           ', num2str(m4)]);
disp(['Skewness Quality per $100, (q>0)         ', num2str(m5)]);
disp(['Kurtosis Quality per $100, (q>0)         ', num2str(m6)]);
disp('---------------------------------------------------');
disp(['Fraction consumed drugs during last 30d  ', num2str(f30d)]);
disp(['Fraction last purchase regular seller    ', num2str(lastreg)]);
disp(['Avg. # of purchases, matched buyers      ', num2str(meanp_m)]);
disp(['Avg. # of purchases, unmatched buyers    ', num2str(meanp_u)]);
disp(['St.dev of purchases, matched buyers      ', num2str(std_m)]);
disp(['St.dev of purchases, unmatched buyers    ', num2str(std_u)]);
disp(['Median # of purchases, match. buyers     ', num2str(med_m)]);
disp(['Median # of purchases, unmatch. buyers   ', num2str(med_u)]);
disp('---------------------------------------------------');
disp(['Fraction consuming drugs in NSDUH (%)    ', num2str(nsduh)]);
disp(['Arrest rate (percent)                    ', num2str(arrest)]);
disp('==================================================');

%% 11) DISPLAY STATISTIC FOR COUNTER FACTUAL SCENARIOS
format short
disp('=======================================================');
disp('Market statistics                         ');
disp(['Rip-offs (%)                             ', num2str(F_0)]);
disp(['Avg. Quality per $100                    ', num2str(m2_s)]);
disp(['Std.dev. Quality per $100                ', num2str(m3_s)]);
disp(['Active Buyers                            ', num2str(B)]);
disp(['Active Sellers                           ', num2str(S/1000)]);
disp(['Fraction matched buyers                  ', num2str(mean(reg))]);
disp(['Avg. # of purchases per month            ', num2str(meanp_t)]);
disp(['Avg. pure grams consumed per month       ', num2str(meanq_t)]);
disp(['Fraction transactions with q = 0         ', num2str(trans_0q)]);
disp('=======================================================');


