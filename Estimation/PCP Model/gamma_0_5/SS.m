% THIS .M FILE EVALUATES STEADY STATE CONDITIONS AT PARAMETER VALUES
% Model 3

function out = SS(x0)

global beta alpha alphatilde sigma sigma_c varphi kappa gamma gamma_H gamma_U...
    gamma_R i_star psi B_bar a eta mu_W thetaP thetaW lambdaW mubar...
    p_star_U p_star_R wd_U wd_R mc_U_U mc_R_R...
    theta_U_UH theta_H_UH theta_U_RH theta_H_RH theta_R_RH...
    theta_U_HU theta_H_HU theta_U_HR theta_H_HR theta_R_HR...
    rho_xi rho_m rho_eps phi_M zeta varepsilon xi_bar a_bar ;



c=x0(1);
n=x0(2);
x=x0(3);
y=x0(4);
Pi=x0(5);
w=x0(6);
mc=x0(7);
q=x0(8);
y_UU=x0(9);
y_UG = x0(10);
y_UR = x0(11);

out(1) = q-0;
out(2) = w-mu_W-log(kappa)-sigma_c*c-varphi*n;
out(3) = mc+log(alphatilde)-(1-alpha)*w-log(1+zeta*i_star)+a_bar;
out(4) = log(1-alpha)+y-n-w+mc;
out(5) = log(alpha)+y-x+mc;
out(6) = exp(y)-exp(y_UU)-exp(y_UG)-exp(y_UR);
out(7) = Pi-(1-exp(mc))*exp(y);
out(8) = log(sigma/(sigma-1))+mc-0;
out(9) = y_UU -log(gamma*(exp(c)+exp(x)));
out(10) = y_UG -log(((1-gamma)/2)*(exp(c)+exp(x)));
out(11) = y_UR -log(((1-gamma)/2)*(exp(c)+exp(x)));
