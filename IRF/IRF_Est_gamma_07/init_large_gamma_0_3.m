clear all;
global beta alpha alphatilde sigma sigma_c varphi kappa gamma ...
     i_star psi B_bar eta mu_W theta_p theta_W lambda_W mubar...
    rho_m rho_eps phi_M zeta varepsilon Gamma mu_bar stdm...
    a_bar rho_a std_a  rho_uip std_uip phi_Y yflex habit;
%params value
beta        = 0.99;                 % Discount factor
sigma       = 2;                    % Demand elasticity 2
sigma_c     = 2;                    % Risk aversion 
varphi      = 2;                    % Frisch elasticity of labor supply 1/varphi =0.5
kappa       = 1;                    % Disutility of labor supply
gamma       = 0.3;                  % home bias
i_star      = 1/beta-1;               % International (gross) interest rate
psi         = 0.001; %0.001;                % interest elasticity of debt 0.001
B_bar       = 0;                    % Steady state debt
eta         = 4;                    % elasticity of substitution across labor varieties
mu_W        = log(eta/(eta-1));     % log SS Wage mark-up
theta_p     = 0.75;                 % price stickiness 0.75
theta_W     = 0.85;                  % Wage rigidity 0.85
lambda_W    = (1-beta*theta_W)*(1-theta_W)/((1+varphi*eta)*theta_W); 
mubar       = log(sigma/(sigma-1)); % log steady state markup
mu_bar      = log(sigma/(sigma-1));     %log of SS markup
zeta        = 0.;                      % interest rate sensitivity of marginal cost
habit       =0.;                  % habit stock as a fraction of lagged consumption H=hC_{-1}; h=0.55



% shocks to TFP
a_bar      = 1.;                    % mean of a 1.9
rho_a      = 0.8; %0.9;                  % persistence 0.71 (what worked, rho_a=0.95, std_a=0.01, correl_a_xi=-0.8
std_a      = 0.01 %0.0539 %0.01 %0.012;                   % SE of eps_a , 0.012


% Monetary Process
rho_m       = 0.5;                     % inertia in interest rates 0.7
rho_eps     = 0.5;                      % persistence in shock to interest rates 0.6
phi_M       = 1.5 %1.5;                      % inflation sensitivity 1.5
phi_Y       = 0.5/4;                        % output sensitivity
stdm        = 0.0025;                    % std. error 

% uip violation
rho_uip =0.5;
std_uip = 0.01;




% pricing regimes

alpha=2/3;
varepsilon=1;
Gamma= varepsilon/(sigma-1);
alphatilde  = alpha^alpha*(1-alpha)^(1-alpha);

%steadystate


x0(1)=0;
x0(2)=log(1-alpha)+x0(1)-0-log(sigma/(sigma-1));
x0(3)=x0(1);
x0(4)=x0(1);
x0(5)=(1-exp(-log(sigma/(sigma-1))))*exp(x0(4));
x0(6)= mu_W+log(kappa)+sigma_c*x0(1)-varphi*x0(2);
x0(7)=-log(sigma/(sigma-1));
x0(8)=0;
x0(9)=log(gamma*(exp(x0(1))+exp(x0(3))));
x0(10)=log(((1-gamma)/2)*(exp(x0(1))+exp(x0(3))));
x0(11)=x0(10);



% Solving for Steady State
options=optimset('MaxFunEvals',1e+8,'MaxIter',10000, 'TolFun', 1e-10);
[x_SS,fval_SS,exitflag]=fsolve(@(x) SS(x), x0,options)
save('init_gamma_0_3.mat', 'x_SS');