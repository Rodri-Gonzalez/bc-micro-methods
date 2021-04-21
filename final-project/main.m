%**************************************************************************
%**************************************************************************
% EMPIRICAL METHODS - FINAL PROJECT 
% Rodrigo Gonzalez
% Simulation of Galenianos & Gavazza (2017)
% "A Structural Model of the Retail Market for Illicit Drugs"
%**************************************************************************
%**************************************************************************
clear all 
cd 'D:\GitDir\bc-micro-methods\final-project'
rng('default');     % For reproductibility

%% GENERATE SAMPLE OF POTENTIAL BUYERS AND SELLERS FROM RECOVERED DISTRIBUTIONS 
param   = parameters()                   ;  % Call parameters
ngrid   = 5000                           ;  % Number of simulations
% BUYERS
mu_z    = param.mu_z                     ;  % Mean valuation
sigz    = param.sigz                     ;  % Stdv of valuation
B_max   = param.B_max                    ;  % Eq number of consumers
p       = param.p                        ;  % Given price
z       = exp(mu_z + sigz*randn(ngrid,1));  % Random draws of consumers tastes
zst     = param.zst                      ;  % Mg utility threshold
R       = p./z                           ;  % Reservation qualities                                  
R_lb    = param.R_lb                     ;  % R lower bar (p/z_ub)
% SELLERS
xi      = param.xi                       ;  % c distribution shape parameter
cst     = param.cst                      ;  % Marginal costs of marginal seller (offer q_lb?))

%% ESTIMATE EQUILIBRIUM THRESHOLDS IMPLIED BY THE MODEL [z*,F_0,R_lb,R_ub,q_lb,q_ub,c_lb,c_ub]
% Set initial values [z* , F_0 , c_ub , q_lb]
% initial =            [148, 0.15, 0.6] ;
% Define function
% thresh  = @(x) th(x, ngrid);
% Search for values that rationalize the parameters
% values  = fminsearch(thresh,initial);

%% USE ESTIMATED THRESHOLDS [z*,F_0,R_lb,R_ub,q_lb,q_ub,c_lb,c_ub] TO GET OPTIMAL QUALITY q*(c)
% Set known threshold values [z* , F_0    , c*  , q_lb]
initial =                    [zst, 0.15862, cst , 0.63];
cgrid   = ([1:1:ngrid]/ngrid).^(1/xi).*cst ;
cgrid   = (cgrid.'); 
[q,c,flag,funct] = optimq(cgrid,initial);  
%plot(c,q);
%histogram(c);
%histogram(q);
q_ub = max(q);

%% SIMULATE MEETINGS

% 




