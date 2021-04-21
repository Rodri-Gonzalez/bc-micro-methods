function [thresholds] = th(init)

% DECLARE PARAMETERS
param   = parameters()   ; % Call parameters function

% FIXED
r       = param.r       ;   % Monthly discount rate
p       = param.p       ;   % Normalized price
K_s     = param.K_s     ;   % Mg. seller monthly flow cost (Levitt&Venkatesh 2000)

% ESTIMATED 
% Frequency rates
alph_b  = param.alph_b  ;   % Rate a buyer meets a seller
gamm    = param.gamm    ;   % Rate of consumption of matched buyer
delt    = param.delt    ;   % Matches' destruction rate

% Consumers' preferences distribution and entry decision
mu_z    = param.mu_z    ;   % Mean of mg. ut distribution
sigz    = param.sigz    ;   % St.dv. of mg. ut distribution
K_b     = param.K_b     ;   % Mg. buyer entry condition
lamb    = param.lamb    ;   % Fraction of pop. with z = 0 taste for crack

% Consumers' reservation quality distribution
mu_R    = param.mu_R    ;   % Mean of reservation quality
sigR    = param.sigR    ;   % St.dv of reservation quality

% Producers' costs distribution and entry
c_st    = param.cst     ;   % Mg. seller marginal cost 
xi      = param.xi      ;   % Shape of mg. cost distribution
c_lbar  = param.c_lbar  ;   % c lower bar

% Selection and measurement error 
mu_et   = param.mu_et   ;   % Selection into ADAM mean  
siget   = param.siget   ;   % Selection into ADAM st.dv.
sigep   = param.sigep   ;   % Quality measurement error st.dv
signu   = param.signu   ;   % Purchases measurement error s.tdv




% Marginal buyer equation
mg_b = alph_b*(zst*integral(x,0,qub) + zst*(integral(f,p/x,qub) - p) - K_b) ; % Equation to find marginal buyer 


% find optimal quality q*(c) supplied by each seller
[c_st,q_st] = ode45(q_prime,c_eq,y0);






end