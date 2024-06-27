
% dcp_05 = load('/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Final Code_June_Binder_Suggestions/Estimation/DCP Model/gamma_0_5/DCP_D/Output/DCP_D_results_rigi_without_phiG.mat','oo_');
% dcp_07 = load('/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Final Code_June_Binder_Suggestions/Estimation/DCP Model/gamma_0_7/DCP_D/Output/DCP_D_results_rigi_without_phiG.mat','oo_');
% 
% lcp_05 = load('/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Final Code_June_Binder_Suggestions/Estimation/LCP Model/gamma_0_5/LCP_D/Output/LCP_D_results_rigi_without_phiG.mat','oo_');
% lcp_07 = load('/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Final Code_June_Binder_Suggestions/Estimation/LCP Model/gamma_0_7/LCP_D/Output/LCP_D_results_rigi_without_phiG.mat','oo_');
% 
% pcp_05 = load('/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Final Code_June_Binder_Suggestions/Estimation/PCP Model/gamma_0_5/PCP_D/Output/PCP_D_results_rigi_without_phiG.mat','oo_');
% pcp_07 = load('/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Final Code_June_Binder_Suggestions/Estimation/PCP Model/gamma_0_7/PCP_D/Output/PCP_D_results_rigi_without_phiG.mat','oo_');

dcp_05 = load('/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Final Code_June_Binder_Suggestions/Estimation/DCP Model/gamma_0_5/DCP_D/Output/DCP_D_results_rigi_with_phiG.mat','oo_');
dcp_07 = load('/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Final Code_June_Binder_Suggestions/Estimation/DCP Model/gamma_0_7/DCP_D/Output/DCP_D_results_rigi_with_phiG.mat','oo_');

lcp_05 = load('/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Final Code_June_Binder_Suggestions/Estimation/LCP Model/gamma_0_5/LCP_D/Output/LCP_D_results_rigi_with_phiG.mat','oo_');
lcp_07 = load('/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Final Code_June_Binder_Suggestions/Estimation/LCP Model/gamma_0_7/LCP_D/Output/LCP_D_results_rigi_with_phiG.mat','oo_');

pcp_05 = load('/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Final Code_June_Binder_Suggestions/Estimation/PCP Model/gamma_0_5/PCP_D/Output/PCP_D_results_rigi_with_phiG.mat','oo_');
pcp_07 = load('/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Final Code_June_Binder_Suggestions/Estimation/PCP Model/gamma_0_7/PCP_D/Output/PCP_D_results_rigi_with_phiG.mat','oo_');

%%

FOntsize_= 20;

figure
tcl = tiledlayout(4,2);

nexttile(tcl)
hh=plot(dcp_05.oo_.prior_density.parameters.phi_pi_U(:,1),dcp_05.oo_.prior_density.parameters.phi_pi_U(:,2),'-k','LineWidth',2);
set(hh,'color',[0.7 0.7 0.7]);
hold on;
plot(dcp_05.oo_.posterior_density.parameters.phi_pi_U(:,1),dcp_05.oo_.posterior_density.parameters.phi_pi_U(:,2),'-k','LineWidth',2);
plot(lcp_05.oo_.posterior_density.parameters.phi_pi_U(:,1),lcp_05.oo_.posterior_density.parameters.phi_pi_U(:,2),'-b','LineWidth',2);
plot(pcp_05.oo_.posterior_density.parameters.phi_pi_U(:,1),pcp_05.oo_.posterior_density.parameters.phi_pi_U(:,2),'-g','LineWidth',2);
plot(dcp_07.oo_.posterior_density.parameters.phi_pi_U(:,1),dcp_07.oo_.posterior_density.parameters.phi_pi_U(:,2),'-k','LineWidth',2,'LineStyle','--');
plot(lcp_07.oo_.posterior_density.parameters.phi_pi_U(:,1),lcp_07.oo_.posterior_density.parameters.phi_pi_U(:,2),'-b','LineWidth',2,'LineStyle','--');
plot(pcp_07.oo_.posterior_density.parameters.phi_pi_U(:,1),pcp_07.oo_.posterior_density.parameters.phi_pi_U(:,2),'-g','LineWidth',2,'LineStyle','--');
h=title('$\phi^U_{\pi}$', 'interpreter', 'latex'); h.FontSize=FOntsize_; grid("minor"); ylim tight; xlim tight;


nexttile(tcl)
hh=plot(dcp_05.oo_.prior_density.parameters.phi_pi_G(:,1),dcp_05.oo_.prior_density.parameters.phi_pi_G(:,2),'-k','LineWidth',2);
set(hh,'color',[0.7 0.7 0.7]);
hold on;
plot(dcp_05.oo_.posterior_density.parameters.phi_pi_G(:,1),dcp_05.oo_.posterior_density.parameters.phi_pi_G(:,2),'-k','LineWidth',2);
plot(lcp_05.oo_.posterior_density.parameters.phi_pi_G(:,1),lcp_05.oo_.posterior_density.parameters.phi_pi_G(:,2),'-b','LineWidth',2);
plot(pcp_05.oo_.posterior_density.parameters.phi_pi_G(:,1),pcp_05.oo_.posterior_density.parameters.phi_pi_G(:,2),'-g','LineWidth',2);
plot(dcp_07.oo_.posterior_density.parameters.phi_pi_G(:,1),dcp_07.oo_.posterior_density.parameters.phi_pi_G(:,2),'-k','LineWidth',2,'LineStyle','--');
plot(lcp_07.oo_.posterior_density.parameters.phi_pi_G(:,1),lcp_07.oo_.posterior_density.parameters.phi_pi_G(:,2),'-b','LineWidth',2,'LineStyle','--');
plot(pcp_07.oo_.posterior_density.parameters.phi_pi_G(:,1),pcp_07.oo_.posterior_density.parameters.phi_pi_G(:,2),'-g','LineWidth',2,'LineStyle','--');
h=title('$\phi^G_{\pi}$', 'interpreter', 'latex'); h.FontSize=FOntsize_; grid("minor"); ylim tight; xlim tight;

nexttile(tcl)
hh=plot(dcp_05.oo_.prior_density.parameters.phi_Y_U(:,1),dcp_05.oo_.prior_density.parameters.phi_Y_U(:,2),'-k','LineWidth',2);
set(hh,'color',[0.7 0.7 0.7]);
hold on;
plot(dcp_05.oo_.posterior_density.parameters.phi_Y_U(:,1),dcp_05.oo_.posterior_density.parameters.phi_Y_U(:,2),'-k','LineWidth',2);
plot(lcp_05.oo_.posterior_density.parameters.phi_Y_U(:,1),lcp_05.oo_.posterior_density.parameters.phi_Y_U(:,2),'-b','LineWidth',2);
plot(pcp_05.oo_.posterior_density.parameters.phi_Y_U(:,1),pcp_05.oo_.posterior_density.parameters.phi_Y_U(:,2),'-g','LineWidth',2);
plot(dcp_07.oo_.posterior_density.parameters.phi_Y_U(:,1),dcp_07.oo_.posterior_density.parameters.phi_Y_U(:,2),'-k','LineWidth',2,'LineStyle','--');
plot(lcp_07.oo_.posterior_density.parameters.phi_Y_U(:,1),lcp_07.oo_.posterior_density.parameters.phi_Y_U(:,2),'-b','LineWidth',2,'LineStyle','--');
plot(pcp_07.oo_.posterior_density.parameters.phi_Y_U(:,1),pcp_07.oo_.posterior_density.parameters.phi_Y_U(:,2),'-g','LineWidth',2,'LineStyle','--');
h=title('$\phi^U_Y$', 'interpreter', 'latex'); h.FontSize=FOntsize_; grid("minor"); ylim tight; xlim tight;

nexttile(tcl)
hh=plot(dcp_05.oo_.prior_density.parameters.phi_Y_G(:,1),dcp_05.oo_.prior_density.parameters.phi_Y_G(:,2),'-k','LineWidth',2);
set(hh,'color',[0.7 0.7 0.7]);
hold on;
plot(dcp_05.oo_.posterior_density.parameters.phi_Y_G(:,1),dcp_05.oo_.posterior_density.parameters.phi_Y_G(:,2),'-k','LineWidth',2);
plot(lcp_05.oo_.posterior_density.parameters.phi_Y_G(:,1),lcp_05.oo_.posterior_density.parameters.phi_Y_G(:,2),'-b','LineWidth',2);
plot(pcp_05.oo_.posterior_density.parameters.phi_Y_G(:,1),pcp_05.oo_.posterior_density.parameters.phi_Y_G(:,2),'-g','LineWidth',2);
plot(dcp_07.oo_.posterior_density.parameters.phi_Y_G(:,1),dcp_07.oo_.posterior_density.parameters.phi_Y_G(:,2),'-k','LineWidth',2,'LineStyle','--');
plot(lcp_07.oo_.posterior_density.parameters.phi_Y_G(:,1),lcp_07.oo_.posterior_density.parameters.phi_Y_G(:,2),'-b','LineWidth',2,'LineStyle','--');
plot(pcp_07.oo_.posterior_density.parameters.phi_Y_G(:,1),pcp_07.oo_.posterior_density.parameters.phi_Y_G(:,2),'-g','LineWidth',2,'LineStyle','--');
h=title('$\phi^G_Y$', 'interpreter', 'latex'); h.FontSize=FOntsize_; grid("minor"); ylim tight; xlim tight;

nexttile(tcl)
hh=plot(dcp_05.oo_.prior_density.parameters.theta_W_U(:,1),dcp_05.oo_.prior_density.parameters.theta_W_U(:,2),'-k','LineWidth',2);
set(hh,'color',[0.7 0.7 0.7]);
hold on;
plot(dcp_05.oo_.posterior_density.parameters.theta_W_U(:,1),dcp_05.oo_.posterior_density.parameters.theta_W_U(:,2),'-k','LineWidth',2);
plot(lcp_05.oo_.posterior_density.parameters.theta_W_U(:,1),lcp_05.oo_.posterior_density.parameters.theta_W_U(:,2),'-b','LineWidth',2);
plot(pcp_05.oo_.posterior_density.parameters.theta_W_U(:,1),pcp_05.oo_.posterior_density.parameters.theta_W_U(:,2),'-g','LineWidth',2);
plot(dcp_07.oo_.posterior_density.parameters.theta_W_U(:,1),dcp_07.oo_.posterior_density.parameters.theta_W_U(:,2),'-k','LineWidth',2,'LineStyle','--');
plot(lcp_07.oo_.posterior_density.parameters.theta_W_U(:,1),lcp_07.oo_.posterior_density.parameters.theta_W_U(:,2),'-b','LineWidth',2,'LineStyle','--');
plot(pcp_07.oo_.posterior_density.parameters.theta_W_U(:,1),pcp_07.oo_.posterior_density.parameters.theta_W_U(:,2),'-g','LineWidth',2,'LineStyle','--');
h=title('$\delta_W^U$', 'interpreter', 'latex'); h.FontSize=FOntsize_; grid("minor"); ylim tight; xlim tight;

nexttile(tcl)
hh=plot(dcp_05.oo_.prior_density.parameters.theta_W_G(:,1),dcp_05.oo_.prior_density.parameters.theta_W_G(:,2),'-k','LineWidth',2);
set(hh,'color',[0.7 0.7 0.7]);
hold on;
plot(dcp_05.oo_.posterior_density.parameters.theta_W_G(:,1),dcp_05.oo_.posterior_density.parameters.theta_W_G(:,2),'-k','LineWidth',2);
plot(lcp_05.oo_.posterior_density.parameters.theta_W_G(:,1),lcp_05.oo_.posterior_density.parameters.theta_W_G(:,2),'-b','LineWidth',2);
plot(pcp_05.oo_.posterior_density.parameters.theta_W_G(:,1),pcp_05.oo_.posterior_density.parameters.theta_W_G(:,2),'-g','LineWidth',2);
plot(dcp_07.oo_.posterior_density.parameters.theta_W_G(:,1),dcp_07.oo_.posterior_density.parameters.theta_W_G(:,2),'-k','LineWidth',2,'LineStyle','--');
plot(lcp_07.oo_.posterior_density.parameters.theta_W_G(:,1),lcp_07.oo_.posterior_density.parameters.theta_W_G(:,2),'-b','LineWidth',2,'LineStyle','--');
plot(pcp_07.oo_.posterior_density.parameters.theta_W_G(:,1),pcp_07.oo_.posterior_density.parameters.theta_W_G(:,2),'-g','LineWidth',2,'LineStyle','--');
h=title('$\delta_W^G$', 'interpreter', 'latex'); h.FontSize=FOntsize_; grid("minor"); ylim tight; xlim tight;


nexttile(tcl)
hh=plot(dcp_05.oo_.prior_density.parameters.theta_p_U(:,1),dcp_05.oo_.prior_density.parameters.theta_p_U(:,2),'-k','LineWidth',2);
set(hh,'color',[0.7 0.7 0.7]);
hold on;
plot(dcp_05.oo_.posterior_density.parameters.theta_p_U(:,1),dcp_05.oo_.posterior_density.parameters.theta_p_U(:,2),'-k','LineWidth',2);
plot(lcp_05.oo_.posterior_density.parameters.theta_p_U(:,1),lcp_05.oo_.posterior_density.parameters.theta_p_U(:,2),'-b','LineWidth',2);
plot(pcp_05.oo_.posterior_density.parameters.theta_p_U(:,1),pcp_05.oo_.posterior_density.parameters.theta_p_U(:,2),'-g','LineWidth',2);
plot(dcp_07.oo_.posterior_density.parameters.theta_p_U(:,1),dcp_07.oo_.posterior_density.parameters.theta_p_U(:,2),'-k','LineWidth',2,'LineStyle','--');
plot(lcp_07.oo_.posterior_density.parameters.theta_p_U(:,1),lcp_07.oo_.posterior_density.parameters.theta_p_U(:,2),'-b','LineWidth',2,'LineStyle','--');
plot(pcp_07.oo_.posterior_density.parameters.theta_p_U(:,1),pcp_07.oo_.posterior_density.parameters.theta_p_U(:,2),'-g','LineWidth',2,'LineStyle','--');
h=title('$\delta_P^U$', 'interpreter', 'latex'); h.FontSize=FOntsize_; grid("minor"); ylim tight; xlim tight;

nexttile(tcl)
hh=plot(dcp_05.oo_.prior_density.parameters.theta_p_G(:,1),dcp_05.oo_.prior_density.parameters.theta_p_G(:,2),'-k','LineWidth',2);
set(hh,'color',[0.7 0.7 0.7]);
hold on;
plot(dcp_05.oo_.posterior_density.parameters.theta_p_G(:,1),dcp_05.oo_.posterior_density.parameters.theta_p_G(:,2),'-k','LineWidth',2);
plot(lcp_05.oo_.posterior_density.parameters.theta_p_G(:,1),lcp_05.oo_.posterior_density.parameters.theta_p_G(:,2),'-b','LineWidth',2);
plot(pcp_05.oo_.posterior_density.parameters.theta_p_G(:,1),pcp_05.oo_.posterior_density.parameters.theta_p_G(:,2),'-g','LineWidth',2);
plot(dcp_07.oo_.posterior_density.parameters.theta_p_G(:,1),dcp_07.oo_.posterior_density.parameters.theta_p_G(:,2),'-k','LineWidth',2,'LineStyle','--');
plot(lcp_07.oo_.posterior_density.parameters.theta_p_G(:,1),lcp_07.oo_.posterior_density.parameters.theta_p_G(:,2),'-b','LineWidth',2,'LineStyle','--');
plot(pcp_07.oo_.posterior_density.parameters.theta_p_G(:,1),pcp_07.oo_.posterior_density.parameters.theta_p_G(:,2),'-g','LineWidth',2,'LineStyle','--');
h=title('$\delta_P^U$', 'interpreter', 'latex'); h.FontSize=FOntsize_; grid("minor"); ylim tight; xlim tight;

% hL = legend('Prior','DCP $\gamma = 0.5$','LCP $\gamma = 0.5$','PCP $\gamma = 0.5$','DCP $\gamma = 0.7$','LCP $\gamma = 0.7$','PCP $\gamma = 0.7$','Orientation','horizontal', 'interpreter', 'latex');
hL = legend({'','DCP $\gamma = 0.5$','LCP $\gamma = 0.5$','PCP $\gamma = 0.5$','DCP $\gamma = 0.7$','LCP $\gamma = 0.7$','PCP $\gamma = 0.7$'},...
    'Orientation','horizontal', 'interpreter', 'latex','NumColumns',3);hL.Layout.Tile = 'South';
set(hL,'FontSize',20);
legend('boxoff')


%% 
% Get the HPD values
% For instance 90% HDP of phi_pi_U for DCP with gamma=0.5 is given by
% dcp_05.oo_.posterior_hpdinf.parameters.phi_pi_U & dcp_05.oo_.posterior_hpdsup.parameters.phi_pi_U

dcp_05.hpd = [dcp_05.oo_.posterior_hpdinf.parameters(1)  dcp_05.oo_.posterior_hpdsup.parameters(1)];
dcp_07.hpd = [dcp_07.oo_.posterior_hpdinf.parameters(1)  dcp_07.oo_.posterior_hpdsup.parameters(1)];

pcp_05.hpd = [pcp_05.oo_.posterior_hpdinf.parameters(1)  pcp_05.oo_.posterior_hpdsup.parameters(1)];
pcp_07.hpd = [pcp_07.oo_.posterior_hpdinf.parameters(1)  pcp_07.oo_.posterior_hpdsup.parameters(1)];

lcp_05.hpd = [lcp_05.oo_.posterior_hpdinf.parameters(1)  lcp_05.oo_.posterior_hpdsup.parameters(1)];
lcp_07.hpd = [lcp_07.oo_.posterior_hpdinf.parameters(1)  lcp_07.oo_.posterior_hpdsup.parameters(1)];


% shocks
dcp_05.hpd_se = [dcp_05.oo_.posterior_hpdinf.shocks_std(1)  dcp_05.oo_.posterior_hpdsup.shocks_std(1)];
dcp_07.hpd_se = [dcp_07.oo_.posterior_hpdinf.shocks_std(1)  dcp_07.oo_.posterior_hpdsup.shocks_std(1)];

pcp_05.hpd_se = [pcp_05.oo_.posterior_hpdinf.shocks_std(1)  pcp_05.oo_.posterior_hpdsup.shocks_std(1)];
pcp_07.hpd_se = [pcp_07.oo_.posterior_hpdinf.shocks_std(1)  pcp_07.oo_.posterior_hpdsup.shocks_std(1)];

lcp_05.hpd_se = [lcp_05.oo_.posterior_hpdinf.shocks_std(1)  lcp_05.oo_.posterior_hpdsup.shocks_std(1)];
lcp_07.hpd_se = [lcp_07.oo_.posterior_hpdinf.shocks_std(1)  lcp_07.oo_.posterior_hpdsup.shocks_std(1)];
