
%--------------------------------------------------------------------------
% plotdata.m
% 
% This file plots the data used for estimation
%
% Necessary files: Model_L_data_load.m
%                  
%
% Last modified: 13.02.2022 (SG)
%--------------------------------------------------------------------------

clear
clc

%% Data for Estimation

load Dataset220324_erdiv100.mat

l=length(FFR(1:75));

DATES = zeros(l,1);

for i = 0:75
   DATES(i+1) = 1999 + i*0.25; 
end
 DATES(76) = [];   
      
% Plot

    figure()
    tiledlayout(3,3) % make it full screen
    nexttile;
    plot(DATES,FFR, 'k','LineWidth',2), title('Federal Funds Rate - U.S.'),         xlim([DATES(1),DATES(end)]), grid on, grid minor
    nexttile;
    plot(DATES,GDP_Def_US,     'k','LineWidth',2), title('Inflation Rate - U.S.'),     xlim([DATES(1),DATES(end)]), grid on, grid minor
    nexttile;
    plot(DATES,GDP_US_Cycle,     'k','LineWidth',2), title('Real GDP Growth - U.S.'), xlim([DATES(1),DATES(end)]), grid on, grid minor
    nexttile;
    plot(DATES,Euri,     'k','LineWidth',2), title('Nominal Interest Rate - Euro Area'),              xlim([DATES(1),DATES(end)]), grid on, grid minor
    nexttile;
    plot(DATES,GDP_Def,      'k','LineWidth',2), title('Inflation Rate - Euro Area'),               xlim([DATES(1),DATES(end)]), grid on, grid minor
    nexttile;
    plot(DATES,GDP_Cycle,       'k','LineWidth',2), title('Real GDP Growth - Euro Area'),                  xlim([DATES(1),DATES(end)]), grid on, grid minor
    nexttile;
    plot(DATES,Exchange_r,    'k','LineWidth',2), title('Exchange Rate'),             xlim([DATES(1),DATES(end)]), grid on, grid minor
    
    
    
    
    
    
    