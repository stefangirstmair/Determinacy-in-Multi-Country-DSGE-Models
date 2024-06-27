
clear
clc

cd gamma07/
dynare DCP_D;
irfModel_1 = oo_.irfs;
clear oo_

dynare PCP_D;
irfModel_2 = oo_.irfs;
clear oo_

dynare LCP_D;
irfModel_3 = oo_.irfs;
clear oo_

cd ..


lengthH = 20;
z = 1:lengthH;
% close all

%%
%%
%%        First Plot
%%
%%


%%  Plotting IRF 

fig(1)= figure('Name',['Monetary Policy Shock in U'],'NumberTitle','off');
% Add a title for the entire figure
sgtitle({'MP Shock in U';sprintf('gamma = %.2f', gamma)});

subplot(3,3,1)
plot(z,irfModel_1.i_G_eps_i_U*100,'LineWidth',2)
hold on
plot(z,irfModel_2.i_G_eps_i_U*100,'LineWidth',2)
plot(z,irfModel_3.i_G_eps_i_U*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Interest Rate in G')
legend('DCP','PCP','LCP')

subplot(3,3,2)
plot(z,irfModel_1.pi_G_eps_i_U*100,'LineWidth',2)
hold on
plot(z,irfModel_2.pi_G_eps_i_U*100,'LineWidth',2)
plot(z,irfModel_3.pi_G_eps_i_U*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Inflation Rate in G')

subplot(3,3,3)
plot(z,irfModel_1.q_G_eps_i_U*100,'LineWidth',2)
hold on
plot(z,irfModel_2.q_G_eps_i_U*100,'LineWidth',2)
plot(z,irfModel_3.q_G_eps_i_U*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Exchange Rate G to U')

subplot(3,3,4)
plot(z,irfModel_1.y_G_eps_i_U*100,'LineWidth',2)
hold on
plot(z,irfModel_2.y_G_eps_i_U*100,'LineWidth',2)
plot(z,irfModel_3.y_G_eps_i_U*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Output in G')

irf_dp_tot_GU_1 = -irfModel_1.p_GU_eps_i_U + irfModel_1.p_UG_eps_i_U;
irf_dp_tot_GU_2 = -irfModel_2.p_GU_eps_i_U + irfModel_2.p_UG_eps_i_U;
irf_dp_tot_GU_3 = -irfModel_3.p_GU_eps_i_U + irfModel_3.p_UG_eps_i_U;
subplot(3,3,5)
plot(z,irf_dp_tot_GU_1*100,'LineWidth',2)
hold on
plot(z,irf_dp_tot_GU_2*100,'LineWidth',2)
plot(z,irf_dp_tot_GU_3*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('ToT in G to U')

%Import Prices
irf_import_UG_ns_eps_i_U_1 = irfModel_1.p_UG_eps_i_U + irfModel_1.q_G_eps_i_U;
irf_import_UG_ns_eps_i_U_2 = irfModel_2.p_UG_eps_i_U + irfModel_2.q_G_eps_i_U;
irf_import_UG_ns_eps_i_U_3 = irfModel_3.p_UG_eps_i_U + irfModel_3.q_G_eps_i_U;
subplot(3,3,6)
plot(z,irf_import_UG_ns_eps_i_U_1*100,'LineWidth',2)
hold on
plot(z,irf_import_UG_ns_eps_i_U_2*100,'LineWidth',2)
plot(z,irf_import_UG_ns_eps_i_U_3*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Import Prices in G from U')

%Export Prices
irf_export_GU_ns_eps_i_U_1 = irfModel_1.p_GU_eps_i_U + irfModel_1.q_G_eps_i_U;
irf_export_GU_ns_eps_i_U_2 = irfModel_2.p_GU_eps_i_U + irfModel_2.q_G_eps_i_U;
irf_export_GU_ns_eps_i_U_3 = irfModel_3.p_GU_eps_i_U + irfModel_3.q_G_eps_i_U;
subplot(3,3,7)
plot(z,irf_export_GU_ns_eps_i_U_1*100,'LineWidth',2)
hold on
plot(z,irf_export_GU_ns_eps_i_U_2*100,'LineWidth',2)
plot(z,irf_export_GU_ns_eps_i_U_3*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Export Prices from G to U')

import_G_1=(irfModel_1.y_RG_eps_i_U+irfModel_1.y_UG_eps_i_U)/2;
import_G_2=(irfModel_2.y_RG_eps_i_U+irfModel_2.y_UG_eps_i_U)/2;
import_G_3=(irfModel_3.y_RG_eps_i_U+irfModel_3.y_UG_eps_i_U)/2;
subplot(3,3,8)
plot(z,import_G_1*100,'LineWidth',2)
hold on
plot(z,import_G_2*100,'LineWidth',2)
plot(z,import_G_3*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Imports G')

export_G_1=(irfModel_1.y_GR_eps_i_U+irfModel_1.y_RU_eps_i_U)/2;
export_G_2=(irfModel_2.y_GR_eps_i_U+irfModel_2.y_RU_eps_i_U)/2;
export_G_3=(irfModel_3.y_GR_eps_i_U+irfModel_3.y_RU_eps_i_U)/2;
subplot(3,3,9)
plot(z,export_G_1*100,'LineWidth',2)
hold on
plot(z,export_G_2*100,'LineWidth',2)
plot(z,export_G_3*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Exports G')


%%  Plotting IRF 

fig(2)= figure('Name',['Monetary Policy Shock in G'],'NumberTitle','off');
sgtitle({'MP Shock in G';sprintf('gamma = %.2f', gamma)});
subplot(3,3,1)
plot(z,irfModel_1.i_G_eps_i_G*100,'LineWidth',2)
hold on
plot(z,irfModel_2.i_G_eps_i_G*100,'LineWidth',2)
plot(z,irfModel_3.i_G_eps_i_G*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Interest Rate in G')
legend('DCP','PCP','LCP')

subplot(3,3,2)
plot(z,irfModel_1.pi_G_eps_i_G*100,'LineWidth',2)
hold on
plot(z,irfModel_2.pi_G_eps_i_G*100,'LineWidth',2)
plot(z,irfModel_3.pi_G_eps_i_G*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Inflation Rate in G')

subplot(3,3,3)
plot(z,irfModel_1.q_G_eps_i_G*100,'LineWidth',2)
hold on
plot(z,irfModel_2.q_G_eps_i_G*100,'LineWidth',2)
plot(z,irfModel_3.q_G_eps_i_G*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Exchange Rate G to U')

subplot(3,3,4)
plot(z,irfModel_1.y_G_eps_i_G*100,'LineWidth',2)
hold on
plot(z,irfModel_2.y_G_eps_i_G*100,'LineWidth',2)
plot(z,irfModel_3.y_G_eps_i_G*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Output in G')

irf_dp_tot_GU_1 = -irfModel_1.p_GU_eps_i_G + irfModel_1.p_UG_eps_i_G;
irf_dp_tot_GU_2 = -irfModel_2.p_GU_eps_i_G + irfModel_2.p_UG_eps_i_G;
irf_dp_tot_GU_3 = -irfModel_3.p_GU_eps_i_G + irfModel_3.p_UG_eps_i_G;
subplot(3,3,5)
plot(z,irf_dp_tot_GU_1*100,'LineWidth',2)
hold on
plot(z,irf_dp_tot_GU_2*100,'LineWidth',2)
plot(z,irf_dp_tot_GU_3*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('ToT in G to U')

%Import Prices
irf_import_UG_ns_eps_i_G_1 = irfModel_1.p_UG_eps_i_G + irfModel_1.q_G_eps_i_G;
irf_import_UG_ns_eps_i_G_2 = irfModel_2.p_UG_eps_i_G + irfModel_2.q_G_eps_i_G;
irf_import_UG_ns_eps_i_G_3 = irfModel_3.p_UG_eps_i_G + irfModel_3.q_G_eps_i_G;
subplot(3,3,6)
plot(z,irf_import_UG_ns_eps_i_G_1*100,'LineWidth',2)
hold on
plot(z,irf_import_UG_ns_eps_i_G_2*100,'LineWidth',2)
plot(z,irf_import_UG_ns_eps_i_G_3*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Import Prices in G from U')

%Export Prices
irf_export_GU_ns_eps_i_G_1 = irfModel_1.p_GU_eps_i_G + irfModel_1.q_G_eps_i_G;
irf_export_GU_ns_eps_i_G_2 = irfModel_2.p_GU_eps_i_G + irfModel_2.q_G_eps_i_G;
irf_export_GU_ns_eps_i_G_3 = irfModel_3.p_GU_eps_i_G + irfModel_3.q_G_eps_i_G;
subplot(3,3,7)
plot(z,irf_export_GU_ns_eps_i_G_1*100,'LineWidth',2)
hold on
plot(z,irf_export_GU_ns_eps_i_G_2*100,'LineWidth',2)
plot(z,irf_export_GU_ns_eps_i_G_3*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Export Prices from G to U')

import_G_1=(irfModel_1.y_RG_eps_i_G+irfModel_1.y_UG_eps_i_G)/2;
import_G_2=(irfModel_2.y_RG_eps_i_G+irfModel_2.y_UG_eps_i_G)/2;
import_G_3=(irfModel_3.y_RG_eps_i_G+irfModel_3.y_UG_eps_i_G)/2;
subplot(3,3,8)
plot(z,import_G_1*100,'LineWidth',2)
hold on
plot(z,import_G_2*100,'LineWidth',2)
plot(z,import_G_3*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Imports G')

export_G_1=(irfModel_1.y_GR_eps_i_G+irfModel_1.y_RU_eps_i_G)/2;
export_G_2=(irfModel_2.y_GR_eps_i_G+irfModel_2.y_RU_eps_i_G)/2;
export_G_3=(irfModel_3.y_GR_eps_i_G+irfModel_3.y_RU_eps_i_G)/2;
subplot(3,3,9)
plot(z,export_G_1*100,'LineWidth',2)
hold on
plot(z,export_G_2*100,'LineWidth',2)
plot(z,export_G_3*100,'LineWidth',2)
hold off
xlim([1 options_.irf]);
hline = refline(0, 0);  hline.Color='k';
title('Exports G')