function param = parameters()

% FIX PARAMETERS 
param.r         = 0.01          ;   % Monthly discount rate
param.p         = 100           ;   % Normalized price
param.K_s       = 1500          ;   % Mg. seller monthly flow cost (Levitt&Venkatesh 2000)
param.B_bar     = 227062.163    ; % Potential buyers - true pop. over 15yr in 2002 (https://www.census.gov/data/tables/time-series/demo/popest/intercensal-2000-2010-national.html) 

% ESTIMATED PARAMETERS

% Frequency rates
param.alph      = 1.27      ;   % Rate a buyer meets a seller
param.gamm      = 19.4      ;   % Rate of consumption of matched buyer
param.delt      = 0.73      ;   % Matches' destruction rate

% Consumers' preferences distribution, entry decision and estimated
param.mu_z      = 5.118     ;   % Mean of mg. ut distribution
param.sigz      = 0.114     ;   % St.dv. of mg. ut distribution
param.z_ub      = exp(param.mu_z+3*param.sigz);  % Consumer with highest valuation   
param.zst       = 150.71    ;   % Mg utility threshold
param.K_b       = 152.64    ;   % Mg. buyer entry condition
param.lamb      = 0.982     ;   % Fraction of pop. with z = 0 taste for crack
param.B_max     = (1-param.lamb)*param.B_bar  ; % Measure of maximum number of active buyers 

% Consumers' reservation quality distribution
param.mu_R      = log(param.p)-param.mu_z ;  % Mean of reservation quality
param.sigR      = param.sigz              ;  % St.dv of reservation quality
param.R_lb      = param.p/param.z_ub      ;

% Producers' costs distribution and entry
param.cst       = 124.368   ;   % Mg. seller marginal cost 
param.xi        = 20.443    ;   % Shape of mg. cost distribution
param.c_lb      = 0         ;   % c lower bar

% Selection and measurement error 
param.mu_et     = -5.237    ;   % Selection into ADAM mean  
param.siget     = 2.803     ;   % Selection into ADAM st.dv.
param.sigep     = 0.526     ;   % Quality measurement error st.dv
param.signu     = 0.522     ;   % Purchases measurement error s.tdv

