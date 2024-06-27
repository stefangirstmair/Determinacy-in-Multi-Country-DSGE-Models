var 
p_UU, p_GG, p_RR,p_UG, p_UR, p_GU, p_GR, p_RG, p_RU,
  p_UU_reset,  p_GG_reset,p_RR_reset,
    c_U,c_R, c_G, 
    pi_U, pi_R, pi_G,
    i_U, i_R, i_G,
    i_U_G, i_U_R,
    B_U, B_G,B_R,
    w_U, w_G, w_R,
    n_U, n_G, n_R,
    mc_U, mc_G, mc_R,
    x_U, x_G, x_R,
    y_U, y_R, y_G,
    y_UU, y_UG, y_UR, y_GG, y_GU, y_GR, y_RR, y_RU, y_RG,
    Profit_U, Profit_G, Profit_R,
    q_G, q_R, 
    eps_M_U, eps_M_G, eps_M_R, a_U, a_G, a_R, uip_G, uip_R
    y_U_UG, y_G_UG, p_U_UG, p_G_UG, p_U_UG_reset, p_G_UG_reset,
    y_U_UR, y_R_UR, p_U_UR, p_R_UR, p_U_UR_reset, p_R_UR_reset,
    y_U_GU, y_G_GU, p_U_GU, p_G_GU, p_U_GU_reset, p_G_GU_reset,
    y_R_GR, y_G_GR, y_U_GR, p_R_GR, p_G_GR, p_U_GR, p_R_GR_reset, p_G_GR_reset, p_U_GR_reset,
    y_U_RU, y_R_RU, p_U_RU, p_R_RU, p_U_RU_reset, p_R_RU_reset,
    y_R_RG, y_G_RG, y_U_RG, p_R_RG, p_G_RG, p_U_RG, p_R_RG_reset, p_G_RG_reset, p_U_RG_reset


    %observables
    d_y_U, d_i_U, d_pi_U, d_w_U
    d_y_G, d_i_G, d_pi_G, d_w_G
    d_q_G
   
    ;
    
    % 
     % p_UG_reset, p_UR_reset, p_GU_reset, p_GR_reset, p_RG_reset,  p_RU_reset, 

varexo eps_i_U, eps_i_G, eps_i_R, eps_R, eps_G, eps_a_U, eps_a_G, eps_a_R, eps_uip_G, eps_uip_R;

parameters beta alpha alphatilde sigma sigma_c varphi kappa gamma 
    i_star psi B_bar eta mu_W theta_p theta_W lambda_W  
     rho_m rho_eps phi_M zeta varepsilon Gamma mu_bar stdm a_bar rho_a
    std_a rho_uip std_uip  phi_Y yflex habit
    
    
    
theta_U_UG
theta_G_UG
theta_U_UR
theta_R_UR
theta_U_GU
theta_G_GU
theta_R_GR
theta_G_GR
theta_U_GR
theta_U_RU
theta_R_RU
theta_R_RG
theta_G_RG
theta_U_RG

phi_pi_U
phi_pi_G
phi_pi_R

phi_Y_U
phi_Y_G
phi_Y_R

rho_m_U
rho_m_G

rho_eps_U
rho_eps_G
rho_eps_R
rho_uip_G
rho_uip_R

gamma_U
gamma_G
    ;
load init_gamma_0_7.mat % Load steady state
% //addpath('init_large.m');
beta        = 0.99;                 % Discount factor
sigma       = 2;                    % Demand elasticity 2
sigma_c     = 2;                    % Risk aversion 
varphi      = 2;                    % Frisch elasticity of labor supply 1/varphi =0.5
kappa       = 1;                    % Disutility of labor supply
gamma       = 0.7;                  % home bias
i_star      = 1/beta-1;               % International (gross) interest rate
psi         = 0.001;%0.000001;  %0.; %               % interest elasticity of debt 0.001
B_bar       = 0;                    % Steady state debt
eta         = 4;                    % elasticity of substitution across labor varieties
mu_W        = log(eta/(eta-1));     % log SS Wage mark-up
theta_p     = 0.75;                 % price stickiness 0.75
theta_W     = 0.85;                  % Wage rigidity 0.85
lambda_W    = (1-beta*theta_W)*(1-theta_W)/((1+varphi*eta)*theta_W); 
mubar       = log(sigma/(sigma-1)); % log steady state markup
mu_bar      = log(sigma/(sigma-1));     %log of SS markup
zeta        = 0;                      % interest rate sensitivity of marginal cost
habit       =0.;                  % habit stock as a fraction of lagged consumption H=hC_{-1}; h=0.55



% shocks to TFP
a_bar      = 1.;                    % mean of a 1.9
rho_a      = 0.8; %0.9;                  % persistence 0.71 (what worked, rho_a=0.95, std_a=0.01, correl_a_xi=-0.8
std_a      = 0.01; %0.0539 %0.01 %0.012;                   % SE of eps_a , 0.012


% Monetary Process
rho_m       = 0.5;%0.5;                     % inertia in interest rates 0.7
rho_eps     = 0.5;                      % persistence in shock to interest rates 0.6
phi_M       = 1.5; %1.5;                      % inflation sensitivity 1.5
phi_Y       = 0.5/4; %0.2;%         %0.00;%              % output sensitivity
stdm        = 0.0025;                    % std. error 

% uip violation
rho_uip =0.5;
std_uip = 0.01;




% pricing regimes

alpha=2/3;
varepsilon=1;
Gamma= varepsilon/(sigma-1);
alphatilde  = alpha^alpha*(1-alpha)^(1-alpha);
yflex=x_SS(4);

% load param1.mat;
% 
% set_param_value('beta',beta);
% set_param_value('alpha',alpha);
% set_param_value('alphatilde',alphatilde);
% set_param_value('sigma',sigma);
% set_param_value('sigma_c',sigma_c);
% set_param_value('varphi',varphi);
% set_param_value('kappa',kappa);
% set_param_value('gamma',gamma);
% set_param_value('i_star',i_star);
% set_param_value('psi',psi);
% set_param_value('B_bar',B_bar);
% set_param_value('eta',eta);
% set_param_value('mu_W',mu_W);
% set_param_value('theta_W',theta_W);
% set_param_value('theta_p',theta_p);
% set_param_value('lambda_W',lambda_W);
% set_param_value('mu_bar',mubar);
% set_param_value('rho_m',rho_m);
% set_param_value('rho_eps',rho_eps);
% set_param_value('phi_M',phi_M);
% set_param_value('zeta',zeta);
% set_param_value('varepsilon',varepsilon);
% set_param_value('Gamma',Gamma);
% set_param_value('stdm',stdm);
% set_param_value('a_bar',a_bar);
% set_param_value('rho_a',rho_a);
% set_param_value('std_a',std_a);
% set_param_value('rho_uip',rho_uip);
% set_param_value('std_uip',std_uip);
% set_param_value('phi_Y',phi_Y);
% set_param_value('yflex',yflex);
% set_param_value('habit',habit);

% % DCP
% theta_U_UG = 1;
% theta_G_UG = 0;
% theta_U_UR = 1;
% theta_R_UR = 0;
% theta_U_GU = 1;
% theta_G_GU = 0;
% theta_R_GR = 0;
% theta_G_GR = 0;
% theta_U_GR = 1;
% theta_U_RU = 1;
% theta_R_RU = 0;
% theta_R_RG = 0;
% theta_G_RG = 0;
% theta_U_RG = 1;

% PCP
theta_U_UG = 1;
theta_G_UG = 0;
theta_U_UR = 1;
theta_R_UR = 0;
theta_U_GU = 0;
theta_G_GU = 1;
theta_R_GR = 0;
theta_G_GR = 1;
theta_U_GR = 0;
theta_U_RU = 0;
theta_R_RU = 1;
theta_R_RG = 1;
theta_G_RG = 0;
theta_U_RG = 0;

% % LCP
% theta_U_UG = 0;
% theta_G_UG = 1;
% theta_U_UR = 0;
% theta_R_UR = 1;
% theta_U_GU = 1;
% theta_G_GU = 0;
% theta_R_GR = 1;
% theta_G_GR = 0;
% theta_U_GR = 0;
% theta_U_RU = 1;
% theta_R_RU = 0;
% theta_R_RG = 0;
% theta_G_RG = 1;
% theta_U_RG = 0;

% % Try
% theta_U_UG = 1; % Works
% theta_G_UG = 0; % Works
% theta_U_UR = 1; % Works
% theta_R_UR = 0; % Works
% theta_U_GU = 1; % Works
% theta_G_GU = 0; % Works 
% theta_R_GR = 1; % Works 
% theta_G_GR = 0; % Works
% theta_U_GR = 0; % Works
% theta_U_RU = 1; % Works
% theta_R_RU = 0; % Works
% theta_R_RG = 0; % Works
% theta_G_RG = 0; % Works
% theta_U_RG = 1; % Works

% phi_pi_U = phi_M;
% phi_pi_G = phi_M;
% phi_pi_R = phi_M;

phi_Y_U = 0.125;
phi_Y_G = 0.125;
phi_Y_R = 0.125;

phi_pi_U = 1.5;
phi_pi_G = 1.5;
phi_pi_R = 1.5;

rho_m_U = 0.5;
rho_m_G = 0.5;

rho_eps_U = 0;%0.5;%
rho_eps_G = 0;%0.5;
rho_eps_R = 0;%0.5;

rho_uip_G = 0.5;
rho_uip_R = 0.5;

gamma_U       = 0.7;                  % home bias
gamma_G       = 0.7;                  % home bias

% phi_pi_U = 0.96;
% phi_pi_G = 0.96;
% phi_pi_R = 1.5;

% THis runs in DCP
% phi_pi_U = 0;%phi_M;%0.88;%phi_M;
% phi_pi_G = 0;%phi_M;%0.06061;%
% phi_pi_R = 5;%phi_M;

% % This runs in PCP
% phi_pi_U = 0.9;
% phi_pi_G = 0;
% phi_pi_R = 0;

model;
// quantities
exp(y_U) = exp(y_UU) + exp(y_UG) + exp(y_UR);
exp(y_G) = exp(y_GG) + exp(y_GU) + exp(y_GR);
exp(y_R) = exp(y_RR) + exp(y_RG) + exp(y_RU);

exp(y_UU) = gamma*(1-varepsilon*(p_UU))^(sigma/varepsilon)*(exp(c_U)+exp(x_U));

exp(y_UG) = theta_U_UG*exp(y_U_UG) + theta_G_UG*exp(y_G_UG);%exp(y_UG) = (1-gamma)/2*(1-varepsilon*(p_UG+q_G))^(sigma/varepsilon)*(exp(c_G)+exp(x_G));
exp(y_U_UG) = (1-gamma)/2*(1-varepsilon*(p_U_UG+q_G))^(sigma/varepsilon)*(exp(c_G)+exp(x_G));
exp(y_G_UG) = (1-gamma)/2*(1-varepsilon*(p_G_UG))^(sigma/varepsilon)*(exp(c_G)+exp(x_G));
p_U_UG_reset = beta*theta_p*(p_U_UG_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_U - Gamma*q_G+ mu_bar);
p_G_UG_reset = beta*theta_p*(p_G_UG_reset(+1)+pi_G(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_U + q_G + mu_bar);
p_U_UG - p_U_UG(-1)+pi_U = (1-theta_p)*(p_U_UG_reset - p_U_UG(-1)+pi_U);
p_G_UG - p_G_UG(-1)+pi_G = (1-theta_p)*(p_G_UG_reset - p_G_UG(-1)+pi_G);

exp(y_UR) = theta_U_UR*exp(y_U_UR) + theta_R_UR*exp(y_R_UR);%exp(y_UR) = (1-gamma)/2*(1-varepsilon*(p_UR+q_R))^(sigma/varepsilon)*(exp(c_R)+exp(x_R));
exp(y_U_UR) = (1-gamma)/2*(1-varepsilon*(p_U_UR+q_R))^(sigma/varepsilon)*(exp(c_R)+exp(x_R));
exp(y_R_UR) = (1-gamma)/2*(1-varepsilon*(p_R_UR))^(sigma/varepsilon)*(exp(c_R)+exp(x_R));
p_U_UR_reset = beta*theta_p*(p_U_UR_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_U - Gamma*q_R+ mu_bar);
p_R_UR_reset = beta*theta_p*(p_R_UR_reset(+1)+pi_R(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_U + q_R + mu_bar);
p_U_UR - p_U_UR(-1)+pi_U = (1-theta_p)*(p_U_UR_reset - p_U_UR(-1)+pi_U);
p_R_UR - p_R_UR(-1)+pi_R = (1-theta_p)*(p_R_UR_reset - p_R_UR(-1)+pi_R);


exp(y_GG) = gamma*(1-varepsilon*(p_GG))^(sigma/varepsilon)*(exp(c_G)+exp(x_G));

exp(y_GU) = theta_U_GU*exp(y_U_GU) + theta_G_GU*exp(y_G_GU);%exp(y_GU) = (1-gamma)/2*(1-varepsilon*(p_GU))^(sigma/varepsilon)*(exp(c_U)+exp(x_U));
exp(y_U_GU) = (1-gamma)/2*(1-varepsilon*(p_U_GU))^(sigma/varepsilon)*(exp(c_U)+exp(x_U));
exp(y_G_GU) = (1-gamma)/2*(1-varepsilon*(p_G_GU-q_G))^(sigma/varepsilon)*(exp(c_U)+exp(x_U));
p_U_GU_reset = beta*theta_p*(p_U_GU_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_G - q_G+ mu_bar);
p_G_GU_reset = beta*theta_p*(p_G_GU_reset(+1)+pi_G(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_G + Gamma*q_G + mu_bar);
p_U_GU - p_U_GU(-1)+pi_U = (1-theta_p)*(p_U_GU_reset - p_U_GU(-1)+pi_U);
p_G_GU - p_G_GU(-1)+pi_G = (1-theta_p)*(p_G_GU_reset - p_G_GU(-1)+pi_G);

exp(y_GR) = theta_R_GR*exp(y_R_GR) + theta_G_GR*exp(y_G_GR) + theta_U_GR*exp(y_U_GR);
exp(y_R_GR) = (1-gamma)/2*(1-varepsilon*(p_R_GR))^(sigma/varepsilon)*(exp(c_R)+exp(x_R));
exp(y_G_GR) = (1-gamma)/2*(1-varepsilon*(p_G_GR+q_R-q_G))^(sigma/varepsilon)*(exp(c_R)+exp(x_R));
exp(y_U_GR) = (1-gamma)/2*(1-varepsilon*(p_U_GR+q_R))^(sigma/varepsilon)*(exp(c_R)+exp(x_R));
p_R_GR_reset = beta*theta_p*(p_R_GR_reset(+1)+pi_R(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_G + q_R - q_G + mu_bar);
p_G_GR_reset = beta*theta_p*(p_G_GR_reset(+1)+pi_G(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_G - Gamma*q_R + mu_bar);
p_U_GR_reset = beta*theta_p*(p_U_GR_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_G- q_G - Gamma*q_R+ mu_bar);
p_R_GR - p_R_GR(-1)+pi_R = (1-theta_p)*(p_R_GR_reset - p_R_GR(-1)+pi_R);
p_G_GR - p_G_GR(-1)+pi_G = (1-theta_p)*(p_G_GR_reset - p_G_GR(-1)+pi_G);
p_U_GR - p_U_GR(-1)+pi_U = (1-theta_p)*(p_U_GR_reset - p_U_GR(-1)+pi_U);


exp(y_RR) = gamma*(1-varepsilon*(p_RR))^(sigma/varepsilon)*(exp(c_R)+exp(x_R));

exp(y_RU) = theta_U_RU*exp(y_U_RU) + theta_R_RU*exp(y_R_RU);%exp(y_RU) = (1-gamma)/2*(1-varepsilon*(p_RU))^(sigma/varepsilon)*(exp(c_U)+exp(x_U));
exp(y_U_RU) = (1-gamma)/2*(1-varepsilon*(p_U_RU))^(sigma/varepsilon)*(exp(c_U)+exp(x_U));
exp(y_R_RU) = (1-gamma)/2*(1-varepsilon*(p_R_RU-q_R))^(sigma/varepsilon)*(exp(c_U)+exp(x_U));
p_U_RU_reset = beta*theta_p*(p_U_RU_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_R - q_R+ mu_bar);
p_R_RU_reset = beta*theta_p*(p_R_RU_reset(+1)+pi_R(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_R + Gamma*q_R + mu_bar);
p_U_RU - p_U_RU(-1)+pi_U = (1-theta_p)*(p_U_RU_reset - p_U_RU(-1)+pi_U);
p_R_RU - p_R_RU(-1)+pi_R = (1-theta_p)*(p_R_RU_reset - p_R_RU(-1)+pi_R);


exp(y_RG) = theta_R_RG*exp(y_R_RG) + theta_G_RG*exp(y_G_RG) + theta_U_RG*exp(y_U_RG);%exp(y_RG) = (1-gamma)/2*(1-varepsilon*(p_RG+q_G))^(sigma/varepsilon)*(exp(c_G)+exp(x_G));
exp(y_R_RG) = (1-gamma)/2*(1-varepsilon*(p_R_RG-q_R+q_G))^(sigma/varepsilon)*(exp(c_G)+exp(x_G));
exp(y_G_RG) = (1-gamma)/2*(1-varepsilon*(p_G_RG))^(sigma/varepsilon)*(exp(c_G)+exp(x_G));
exp(y_U_RG) = (1-gamma)/2*(1-varepsilon*(p_U_RG+q_G))^(sigma/varepsilon)*(exp(c_G)+exp(x_G));
p_R_RG_reset = beta*theta_p*(p_R_RG_reset(+1)+pi_R(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_R - Gamma*q_G + mu_bar);
p_G_RG_reset = beta*theta_p*(p_G_RG_reset(+1)+pi_G(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_R - q_R + q_G + mu_bar); % Not sure, also works with (mc_R - q_R  + mu_bar)
p_U_RG_reset = beta*theta_p*(p_U_RG_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_R- q_R - Gamma*q_G+ mu_bar);
p_R_RG - p_R_RG(-1)+pi_R = (1-theta_p)*(p_R_RG_reset - p_R_RG(-1)+pi_R);
p_G_RG - p_G_RG(-1)+pi_G = (1-theta_p)*(p_G_RG_reset - p_G_RG(-1)+pi_G);
p_U_RG - p_U_RG(-1)+pi_U = (1-theta_p)*(p_U_RG_reset - p_U_RG(-1)+pi_U);



// CPI 
% gamma*p_UU+((1-gamma)/2)*(p_GU+p_RU)=0;             
% gamma*p_GG+((1-gamma)/2)*(p_UG+q_G+p_RG+q_G)=0;
% gamma*p_RR+((1-gamma)/2)*(p_UR+q_R+p_GR+q_R)=0;
gamma*p_UU+((1-gamma)/2)*(theta_U_GU*(p_U_GU) + theta_G_GU*(p_G_GU-q_G) + theta_U_RU*(p_U_RU)+ theta_R_RU*(p_R_RU-q_R))=0;
gamma*p_GG+((1-gamma)/2)*(theta_U_UG*(p_U_UG+q_G) + theta_G_UG*(p_G_UG) + theta_U_RG*(p_U_RG+q_G) + theta_R_RG*(p_R_RG-q_R+q_G) + theta_G_RG*(p_G_RG))=0;
gamma*p_RR+((1-gamma)/2)*(theta_U_UR*(p_U_UR+q_R) + theta_R_UR*(p_R_UR) + theta_R_GR*(p_R_GR) + theta_G_GR*(p_G_GR+q_R-q_G) + theta_U_GR*(p_U_GR+q_R))=0;



% // Price reset equations
p_UU_reset = beta*theta_p*(p_UU_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_U + mu_bar);
p_GG_reset = beta*theta_p*(p_GG_reset(+1)+pi_G(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_G + mu_bar);
p_RR_reset = beta*theta_p*(p_RR_reset(+1)+pi_R(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_R + mu_bar);
% 
% p_UG_reset = beta*theta_p*(p_UG_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_U - Gamma*q_G+ mu_bar);
% p_UR_reset = beta*theta_p*(p_UR_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_U - Gamma*q_R+ mu_bar);
% 
% p_GU_reset = beta*theta_p*(p_GU_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_G - q_G+ mu_bar);
% p_GR_reset = beta*theta_p*(p_GR_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_G- q_G - Gamma*q_R+ mu_bar);
% 
% p_RU_reset = beta*theta_p*(p_RU_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_R - q_R+ mu_bar);
% p_RG_reset = beta*theta_p*(p_RG_reset(+1)+pi_U(+1)) + (1-beta*theta_p)/(1+Gamma)*(mc_R- q_R - Gamma*q_G+ mu_bar);


% // Price evolution
p_UU - p_UU(-1)+pi_U = (1-theta_p)*(p_UU_reset - p_UU(-1)+pi_U);
p_GG - p_GG(-1)+pi_G = (1-theta_p)*(p_GG_reset - p_GG(-1)+pi_G);
p_RR - p_RR(-1)+pi_R = (1-theta_p)*(p_RR_reset - p_RR(-1)+pi_R);
% 
% p_UG - p_UG(-1)+pi_U = (1-theta_p)*(p_UG_reset - p_UG(-1)+pi_U);
% p_UR - p_UR(-1)+pi_U = (1-theta_p)*(p_UR_reset - p_UR(-1)+pi_U);
% 
% p_GU - p_GU(-1)+pi_U = (1-theta_p)*(p_GU_reset - p_GU(-1)+pi_U);
% p_GR - p_GR(-1)+pi_U = (1-theta_p)*(p_GR_reset - p_GR(-1)+pi_U);
% 
% p_RU - p_RU(-1)+pi_U = (1-theta_p)*(p_RU_reset - p_RU(-1)+pi_U);
% p_RG - p_RG(-1)+pi_U = (1-theta_p)*(p_RG_reset - p_RG(-1)+pi_U);

p_UR = theta_U_UR*(p_U_UR) + theta_R_UR*(p_R_UR);
p_UG = theta_U_UG*(p_U_UG) + theta_G_UG*(p_G_UG);
p_GU = theta_U_GU*(p_U_GU) + theta_G_GU*(p_G_GU);
p_GR = theta_R_GR*(p_R_GR) + theta_G_GR*(p_G_GR) + theta_U_GR*(p_U_GR);
p_RU = theta_U_RU*(p_U_RU) + theta_R_RU*(p_R_RU);
p_RG = theta_R_RG*(p_R_RG) + theta_G_RG*(p_G_RG) + theta_U_RG*(p_U_RG);


// Euler equations
exp(c_U-habit*c_U(-1))^(-sigma_c) = beta*(1+i_U)*exp(c_U(+1)-habit*c_U)^(-sigma_c)/(1+pi_U(+1));
exp(c_G-habit*c_G(-1))^(-sigma_c) = beta*(1+i_G)*exp(c_G(+1)-habit*c_G)^(-sigma_c)/(1+pi_G(+1));
exp(c_R-habit*c_R(-1))^(-sigma_c) = beta*(1+i_R)*exp(c_R(+1)-habit*c_R)^(-sigma_c)/(1+pi_R(+1));

// Budget constraints
exp(c_G) + (exp(q_G)*(1+i_U_G(-1))*B_G(-1)/(1+pi_U)) = Profit_G + exp(q_G)*B_G+exp(w_G+n_G);
exp(c_R) + (exp(q_R)*(1+i_U_R(-1))*B_R(-1)/(1+pi_U)) = Profit_R + exp(q_R)*B_R+exp(w_R+n_R);

// Wages

(w_U-w_U(-1)) + pi_U = beta*(w_U(+1)-w_U + pi_U(+1)) + (1-beta*theta_W)*(1-theta_W)/((1+varphi*eta)*theta_W)*((sigma_c/(1-habit))*(c_U-habit*c_U(-1)) + varphi*n_U + log(kappa) - w_U + mu_W);
(w_G-w_G(-1)) + pi_G = beta*(w_G(+1)-w_G + pi_G(+1)) + (1-beta*theta_W)*(1-theta_W)/((1+varphi*eta)*theta_W)*((sigma_c/(1-habit))*(c_G-habit*c_G(-1)) + varphi*n_G + log(kappa) - w_G + mu_W);
(w_R-w_R(-1)) + pi_R = beta*(w_R(+1)-w_R + pi_R(+1)) + (1-beta*theta_W)*(1-theta_W)/((1+varphi*eta)*theta_W)*((sigma_c/(1-habit))*(c_R-habit*c_R(-1)) + varphi*n_R + log(kappa) - w_R + mu_W);

// Marginal cost

mc_U = (1-alpha)*w_U + log(1+zeta*i_U(-1))- a_U - (alpha*log(alpha)+(1-alpha)*log(1-alpha));
mc_G = (1-alpha)*w_G + log(1+zeta*i_U_G(-1))- a_G - (alpha*log(alpha)+(1-alpha)*log(1-alpha));
mc_R = (1-alpha)*w_R + log(1+zeta*i_U_R(-1))- a_R - (alpha*log(alpha)+(1-alpha)*log(1-alpha));

// labor demand
(1-alpha)*exp(y_U)/exp(n_U) = exp(w_U)/exp(mc_U);
(1-alpha)*exp(y_G)/exp(n_G) = exp(w_G)/exp(mc_G);
(1-alpha)*exp(y_R)/exp(n_R) = exp(w_R)/exp(mc_R);


// Intermediate input demand
alpha*exp(y_U)/exp(x_U) = 1/exp(mc_U);
alpha*exp(y_G)/exp(x_G) = 1/exp(mc_G);
alpha*exp(y_R)/exp(x_R) = 1/exp(mc_R);


% Profit
Profit_U = exp(p_UU + y_UU) + theta_U_UR*exp(p_U_UR + y_U_UR) + theta_R_UR*exp(q_R + p_R_UR + y_R_UR) + theta_U_UG*exp(p_U_UG + y_U_UG) + theta_G_UG*exp(p_G_UG + y_G_UG) - exp(mc_U+y_U);
Profit_G = exp(p_GG + y_GG) + theta_U_GU*exp(q_G + p_U_GU + y_U_GU) + theta_G_GU*exp(p_G_GU + y_G_GU) + theta_U_GR*exp(q_G+p_U_GR+y_U_GR) + theta_G_GR*exp(p_G_GR+y_G_GR) + theta_R_GR*exp(-q_R+p_R_GR+q_G+y_R_GR)-exp(mc_G+y_G);
Profit_R = exp(p_RR + y_RR) + theta_U_RU*exp(q_R + p_U_RU + y_U_RU) + theta_R_RU*exp(p_R_RU + y_R_RU) + theta_R_RG*exp(p_R_RG + y_R_RG) + theta_U_RG*exp(q_R + p_U_RG + y_U_RG) + theta_G_RG*exp(q_R + p_G_RG - q_G + y_G_RG)-exp(mc_R+y_R);

// Dollar bond interest rates in G and R
i_U_G = i_U + psi*(exp(B_G-B_bar)-1) + eps_G;
i_U_R = i_U + psi*(exp(B_R-B_bar)-1) + eps_R;

// Domestic interest rates

i_U-i_star = rho_m_U *(i_U(-1)-i_star) + (1-rho_m_U)*(phi_pi_U*pi_U +phi_Y_U*(y_U-yflex)) - eps_M_U;
i_G-i_star = rho_m_G *(i_G(-1)-i_star) + (1-rho_m_G)*(phi_pi_G*pi_G +phi_Y_G*(y_G-yflex)) - eps_M_G;
i_R-i_star = rho_m *(i_R(-1)-i_star) + (1-rho_m)*(phi_pi_R*pi_R +phi_Y_R*(y_R-yflex)) - eps_M_R;

// exchange rates

i_G-i_U_G = q_G(+1) -q_G + pi_G(+1)-pi_U(+1)+uip_G;
i_R-i_U_R = q_R(+1) -q_R + pi_R(+1)-pi_U(+1)+uip_R;

// bond market clearing
B_U+B_G+B_R=0;

// Shocks to monetary rule
eps_M_U=rho_eps_U*eps_M_U(-1) + eps_i_U;
eps_M_G=rho_eps_G*eps_M_G(-1) + eps_i_G;
eps_M_R=rho_eps_R*eps_M_R(-1) + eps_i_R;

// shocks to UIP
uip_G=rho_uip_G*uip_G(-1)+eps_uip_G;
uip_R=rho_uip_R*uip_R(-1)+eps_uip_R;

// shocks to productivity
a_U-a_bar=rho_a*(a_U(-1)-a_bar)+eps_a_U;
a_G-a_bar=rho_a*(a_G(-1)-a_bar)+eps_a_G;
a_R-a_bar=rho_a*(a_R(-1)-a_bar)+eps_a_R;


// Other variables
//gdp = log(exp(c)+exp(p_HU)*exp(y_HU)+exp(p_HR)*exp(y_HR)+exp(e_U)*exp(xi)-exp(p_UH)*exp(y_UH)-exp(p_RH)*exp(y_RH));
//CAGDP = exp(e_U)*(B_U(-1)-B_U)/exp(gdp);
//TBGDP = (exp(p_HU)*exp(y_HU)+exp(p_HR)*exp(y_HR)+exp(e_U)*exp(xi)-exp(p_UH)*exp(y_UH)-exp(p_RH)*exp(y_RH))/exp(gdp);
//TBMGDP = (exp(p_HU)*exp(y_HU)+exp(p_HR)*exp(y_HR)-exp(p_UH)*exp(y_UH)-exp(p_RH)*exp(y_RH))/exp(gdp);
//labprod=gdp-n;
//manuf_va=log(Profit+exp(w+l));
//exportsM = log(exp(y_HU)+exp(y_HR));
//importsM = log(exp(y_UH)+exp(y_RH));



% Observables
 d_y_U = y_U - steady_state(y_U);
 d_i_U = i_U - steady_state(i_U);
 d_pi_U = pi_U;
 d_w_U = w_U - steady_state(w_U);

 d_y_G = y_G - steady_state(y_G);
 d_i_G = i_G - steady_state(i_G);
 d_pi_G = pi_G;
 d_w_G = w_G - steady_state(w_G);

 d_q_G = q_G;

end;

initval;
c_U=x_SS(1);c_G=x_SS(1); c_R=x_SS(1);
n_U=x_SS(2);n_G=x_SS(2);n_R=x_SS(2);
x_U=x_SS(3);x_G=x_SS(3);x_R=x_SS(3);
y_U=x_SS(4);y_G=x_SS(4);y_R=x_SS(4);
Profit_U=x_SS(5);Profit_G=x_SS(5);Profit_R=x_SS(5);
w_U=x_SS(6);w_G=x_SS(6);w_R=x_SS(6);
mc_U=x_SS(7);mc_G=x_SS(7);mc_R=x_SS(7);
q_G=x_SS(8);q_R=x_SS(8);
y_UU=x_SS(9);y_RR=x_SS(9);y_GG=x_SS(9);
% y_UG=x0(10);y_UR=x0(10);y_GU=x0(10);y_GR=x0(10);y_RU=x0(10);y_RG=x0(10);
y_UG=x_SS(10);y_UR=x_SS(10);y_GU=x_SS(10);y_GR=x_SS(10);y_RU=x_SS(10);y_RG=x_SS(10);
B_G =B_bar;B_U=B_bar;B_R= B_bar;
pi_U=0;pi_G=0;pi_R=0;
i_U =i_star;i_R=i_star;i_G= i_star;
i_U_G =i_star;i_U_R= i_star;

% p_UU_reset=0;p_UG_reset=0;p_UR_reset=0;
% p_RR_reset=0;p_RG_reset=0;p_RU_reset=0;
% p_GG_reset=0;p_GR_reset=0;p_GU_reset=0;

p_UU=0;p_UG=0;p_UR=0;
p_RR=0;p_RG=0;p_RU=0;
p_GG=0;p_GR=0;p_GU=0;

eps_M_U=0; eps_M_R=0; eps_M_G=0;
a_U=a_bar; a_G=a_bar; a_R=a_bar;

y_U_UG=x_SS(10);y_U_UR=x_SS(10);y_U_GU=x_SS(10);y_U_GR=x_SS(10);y_U_RU=x_SS(10);y_U_RG=x_SS(10);
y_G_UG=x_SS(10);y_G_GU=x_SS(10);y_G_GR=x_SS(10);y_G_RG=x_SS(10);
y_R_UR=x_SS(10);y_R_GR=x_SS(10);y_R_RU=x_SS(10);y_R_RG=x_SS(10);

end;
steady(maxit=10000);
% check;
% model_diagnostics;

%  //monetary policy shock in US
 % shocks;
 % var eps_i_U;
 % stderr stdm;
 % end;

% //monetary policy shock in G
% shocks;
% var eps_i_G;
% stderr 0.0025;
% end;

shocks;
var eps_a_G;
stderr 0.01;
end;


//shocks to eU-eR
//shocks;
//var eps_R;
//stderr std_eR;  //0.002
//end;


% stoch_simul(order=1, periods=1000, drop=100, irf=20, nograph) i_G pi_G q_G y_G c_G mc_G; // test for output and interest rate of U % 
% stoch_simul(order=1, periods=1000, drop=100, irf=20) i_G pi_G q_G y_G c_G mc_G; // test for output and interest rate of U % 
    % stoch_simul(order=1, periods=1000, drop=100, irf=20, nograph) pi_G pi_U  pi_R y_G y_U y_R i_G i_U i_R; // test for output and interest rate of U % 

     % stoch_simul(order=1, periods=1000, drop=100, irf=20) pi_G pi_U  pi_R y_G y_U y_R i_G i_U i_R; // test for output and interest rate of U % 
   
% savefile = 'simulation.mat';
% save(savefile, 'oo_'); 



% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% SENSITIVIITY ANALYSIS %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% estimated_params;
% phi_Y_U, 0.125, 0,2;
% phi_Y_G, 0.125, 0,2;
% phi_Y_R, 0.125, 0,2;
% 
% phi_pi_U, 1.5, 0,2;
% phi_pi_G, 1.5, 0,2;
% phi_pi_R, 1.5, 0,2;
% 
% gamma, 0.7, 0, 1;
% end;
% 
% options_.varobs = M_.endo_names; % dynare_sensitivity requires varobs block
%                                  % alternative and quick way to assume all variables are observbable:
%                                  % options_.varobs = M_.endo_names; 
% 
% dynare_sensitivity;


varobs 
d_y_U 
d_i_U 
d_pi_U 
% d_w_U
d_y_G 
d_i_G 
d_pi_G 
% d_w_G
d_q_G 
;
% varobs  d_i_U d_pi_U d_y_U  d_i_G d_pi_G  d_y_G;

% load Dataset011123_3.mat;
% load Dataset010224.mat;
% load Dataset220324.mat;
load Dataset220324_erdiv100.mat;

// TIME SPAN
BEG_OBS          = 1;   
END_OBS          = length(DATES);      
DATES            = DATES(BEG_OBS:END_OBS,:);
datadescr.DATES  = DATES;
datadescr.MODEL  = M_.fname;

dataseries.d_y_U  = GDP_US_Cycle(BEG_OBS:END_OBS);
datadescr.d_y_U   = 'GDP US';

dataseries.d_i_U  = FFR(BEG_OBS:END_OBS);
datadescr.d_i_U   = 'Interest Rate US';

dataseries.d_pi_U  = GDP_Def_US(BEG_OBS:END_OBS);
datadescr.d_pi_U   = 'Inflation Rate US';

dataseries.d_w_U  = Wage_US_Cycle(BEG_OBS:END_OBS);
datadescr.d_w_U   = 'Wages US';

dataseries.d_y_G  = GDP_Cycle(BEG_OBS:END_OBS);
datadescr.d_y_G   = 'GDP EA';

dataseries.d_i_G  = Euri(BEG_OBS:END_OBS);
datadescr.d_i_G   = 'Interest Rate EA';

dataseries.d_pi_G  = GDP_Def(BEG_OBS:END_OBS);
datadescr.d_pi_G   = 'Inflation Rate EA';

dataseries.d_w_G  = Wage_Cycle(BEG_OBS:END_OBS);
datadescr.d_w_G   = 'Wages US';

dataseries.d_q_G  = Exchange_r(BEG_OBS:END_OBS);
datadescr.d_q_G   = 'Exchange Rate';

save ('DYNARE_estimdata','-STRUCT','dataseries')

estimated_params;
  

stderr eps_i_U, 0.0031, 0, inf, inv_gamma_pdf, 0.01, 0.005  ; 
stderr eps_i_G, 0.0025, 0, inf, inv_gamma_pdf, 0.01, 0.005  ; 
stderr eps_i_R, 0.0071, 0, inf, inv_gamma_pdf, 0.01, 0.005  ; 
stderr eps_R, 0.0081, 0, inf, inv_gamma_pdf, 0.01, 0.005  ; 
stderr eps_G, 0.0109, 0, inf, inv_gamma_pdf, 0.01, 0.005  ; 
stderr eps_a_U, 0.0050, 0, inf, inv_gamma_pdf, 0.01, 0.005  ; 
stderr eps_a_G, 0.0051, 0, inf, inv_gamma_pdf, 0.01, 0.005  ; 
stderr eps_a_R, 0.0074, 0, inf, inv_gamma_pdf, 0.01, 0.005  ; 
stderr eps_uip_G, 0.0029, 0, inf, inv_gamma_pdf, 0.01, 0.005  ; 
stderr eps_uip_R, 0.0064, 0, inf, inv_gamma_pdf, 0.01, 0.005  ; 

phi_pi_U , 1.7958, gamma_pdf, 1.5, 0.25;
phi_pi_G , 1.4228, gamma_pdf, 1.5, 0.25;
phi_pi_R , 1.7235, gamma_pdf, 1.5, 0.25;

phi_Y_U , 0.7403, normal_pdf, 0.25, 0.15;
phi_Y_G , 0.7267, normal_pdf, 0.25, 0.15;
phi_Y_R , 0.4536, normal_pdf, 0.25, 0.15;

theta_p, 0.7063, beta_pdf, 0.75, 0.1;
theta_W, 0.6990, beta_pdf, 0.75, 0.1;

rho_m_U, 0.5511, beta_pdf, 0.75, 0.1;
rho_m_G, 0.5585, beta_pdf, 0.75, 0.1;
rho_m, 0.4380, beta_pdf, 0.75, 0.1;

rho_eps_U, 0.9341, beta_pdf, 0.75, 0.1;
rho_eps_G, 0.9101, beta_pdf, 0.75, 0.1;
rho_eps_R, 0.9331, beta_pdf, 0.75, 0.1;
rho_uip_G, 0.9414, beta_pdf, 0.75, 0.1;
rho_uip_R, 0.7279, beta_pdf, 0.75, 0.1;

% gamma_U, , beta_pdf, 0.7, 0.1;
% gamma_G, , beta_pdf, 0.7, 0.1;
end;

% options_.debug=1;

estimated_params_init(use_calibration);
end;
% 
% identification;

estimation(datafile=DYNARE_estimdata,mh_replic=250000, mh_nblocks=1, mh_drop=0.5,plot_priors = 0, mode_compute=4
, prior_trunc=0);  %  , mode_file = DCP_D_Est_mode250324

% estimation(datafile=DYNARE_estimdata,mh_replic=250000, mh_nblocks=1, mh_drop=0.5,plot_priors = 0, mode_compute=4
% , prior_trunc=0, mode_file = DCP_D_Est_mode250324);  
% ESTIMATION RESULTS
% 
% Log data density is 2250.925737.
% 
% parameters
%              prior mean   post. mean        90% HPD interval    prior       pstdev
% 
% phi_pi_U          1.500       1.6076      1.0544      2.1345     gamm       0.2500
% phi_pi_G          1.500       1.3865      0.9529      1.8199     gamm       0.2500
% phi_pi_R          1.500       1.8287      1.3422      2.3188     gamm       0.2500
% phi_Y_U           0.250       0.8221      0.5924      1.0443     norm       0.1500
% phi_Y_G           0.250       0.7789      0.5953      0.9762     norm       0.1500
% phi_Y_R           0.250       0.4717      0.2726      0.6630     norm       0.1500
% theta_p           0.750       0.6938      0.6378      0.7446     beta       0.1000
% theta_W           0.750       0.7352      0.5760      0.9626     beta       0.1000
% rho_m_U           0.750       0.5564      0.4648      0.6473     beta       0.1000
% rho_m_G           0.750       0.5689      0.4741      0.6687     beta       0.1000
% rho_m             0.750       0.4314      0.2914      0.5753     beta       0.1000
% rho_eps_U         0.750       0.9153      0.8727      0.9575     beta       0.1000
% rho_eps_G         0.750       0.8932      0.8464      0.9404     beta       0.1000
% rho_eps_R         0.750       0.9285      0.8981      0.9582     beta       0.1000
% rho_uip_G         0.750       0.9381      0.9230      0.9538     beta       0.1000
% rho_uip_R         0.750       0.7091      0.6095      0.8139     beta       0.1000
% 
% standard deviation of shocks
%              prior mean   post. mean        90% HPD interval    prior       pstdev
% 
% eps_i_U           0.010       0.0031      0.0025      0.0036     invg       0.0050
% eps_i_G           0.010       0.0025      0.0021      0.0029     invg       0.0050
% eps_i_R           0.010       0.0076      0.0050      0.0103     invg       0.0050
% eps_R             0.010       0.0076      0.0050      0.0101     invg       0.0050
% eps_G             0.010       0.0110      0.0089      0.0129     invg       0.0050
% eps_a_U           0.010       0.0056      0.0039      0.0073     invg       0.0050
% eps_a_G           0.010       0.0059      0.0040      0.0078     invg       0.0050
% eps_a_R           0.010       0.0112      0.0045      0.0191     invg       0.0050
% eps_uip_G         0.010       0.0030      0.0025      0.0036     invg       0.0050
% eps_uip_R         0.010       0.0066      0.0048      0.0084     invg       0.0050
% Total computing time : 0h48m28s