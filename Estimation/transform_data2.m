% This file transforms the data such that it can be used in a
% log-linearized DSGE model
% They are transfromed in a way to match model variables one-to-one
% Edited: Stefan Girstmair 16.05.2023
% Chair of International Macroeconomics and Macroeconometrics

% function Zobs = transform_data

% % Housekeeping
% clear
% clc

% Read in data downloaded from AWM Database
% (Data codes are in first line)
% data = xlsread('awm19up18.csv');
data = readmatrix('awm19up18.csv');
data(:,1) = [];

% Read in data downloaded from FRED
% (Data codes are in first line)
% data_US = xlsread('GDP.xls');
data_US = readmatrix('GDP.xls');
data_US(:,1) = [];
% Start date 1970Q1
% Start date 2017Q4

% q = 85;%1990q4
q = 118; %1999q1
tt = 192;

% Beginn at 1983Q1 as Villa (period 53)

% Gross Domestic Product (GDP) at market prices, Million Euro, Chain linked volume, Calendar and seasonally adjusted data, Reference year 1995.
% Transformation:
% Log and HP-Filter to make stationary
% Divided by price level to make real
GDP = log(data(q:tt,1)./data(q:tt,7));
[GDP_Trend,GDP_Cycle] = one_sided_hp_filter(GDP);

% Individual Consumption Expenditure, Millions of euros, Chain linked volume, Calendar and seasonally adjusted data, Reference year 1995.
% Transformation:
% Log and HP-Filter to make stationary
% Divided by price level to make real
Cons = log(data(q:tt,2)./data(q:tt,7));
[Cons_Trend,Cons_Cycle] = one_sided_hp_filter(Cons);

% HICP - Overall Index, Index, Neither seasonally nor working day adjusted data, Index base year 1996 (1996 = 100).
% Transformation:
% Log and First Difference
HICP = log(data(72:tt,19));
HICP = (diff(HICP));
HICP = HICP - mean(HICP);

% GDP Deflator, Index, Index base year 1995 (1995 = 1). Defined as the ratio of nominal, and real gross domestic product (GDP).
% Transformation:
% Log and First Difference
GDP_Def = log(data(q-1:tt,7));
GDP_Def = (diff(GDP_Def));
GDP_Def = GDP_Def - mean(GDP_Def);

% Nominal Short-Term Interest Rate, Euribor 3-month, Percent per annum, Last trade price.
% Transformation:
% Log and Make it Quarterly
Euri = data(q:tt,33);
Euri = log(1+Euri/400) - mean(log(1+Euri/400));

% Gross Fixed Capital Formation, Millions of euros, Chain linked volume, Calendar and seasonally adjusted data, Reference year 1995.
% Transformation:
% Log and HP-Filter to make stationary
Inve = log(data(q:tt,4));
[Inve_Trend,Inve_Cycle] = one_sided_hp_filter(Inve);

% % Wage per Head. Calculated as the ratio of compensation of employees, and total employment (WRN = WIN / LNN).
% % Transformation:
% % Log and HP-Filter to make stationary
% % Divided by price level to make real
Wage = log(data(q:tt,44)./data(q:tt,7));
[Wage_Trend,Wage_Cycle] = one_sided_hp_filter(Wage);
% Wage per Head. Calculated as the ratio of compensation of employees, and total employment (WRN = WIN / LNN).
% Transformation:
% Log and HP-Filter to make stationary
% Divided by price level to make real
% Wage = log(data_US(q:tt,11)./data(q:tt,7));
% [Wage_Trend,Wage_Cycle] = one_sided_hp_filter(Wage);

% Euro-per-USD Exchange Rate.
Exchange = log(data(q:tt,47));
Exchange_n = Exchange - mean(Exchange);
Exchange_n=Exchange_n./100;%./10;

[Exchange_Trend,Exchange_Cycle] = one_sided_hp_filter(Exchange);
Exchange_Cycle=Exchange_Cycle./10;

Exchange_Gr = log(data(q-1:tt,47));
Exchange_Gr = (diff(Exchange_Gr));
Exchange_Gr = Exchange_Gr - mean(Exchange_Gr);
Exchange_Gr=Exchange_Gr./10;

% % Euro-per-USD Exchange Rate. From Fred, spot ER*(CPI/HICP)
Exchange_r = log(1./data_US(q:tt,12)) ;
Exchange_r = (Exchange_r - mean(Exchange_r))./100;

% % Euro-per-USD Exchange Rate. From Fred, spot ER*(CPI/HICP)
Exchange_r_diff = log(1./data_US(q-1:tt,12)) ;
Exchange_r_diff = diff(Exchange_r_diff);

% Euro-Dollar Real Exchange Rate Dynamics in an Estimated Two-Country Model:
% What is Important and What is Not (2005)

% Finally, we construct the real exchange rate series by taking the nominal exchange rate in
% terms of euros per U.S. dollar and multiplying it by the ratio of consumer price indices. The
% “synthetic” euro/U.S. dollar exchange rate prior to the launch of the euro in 1999 also comes
% from Eurostat, while the US CPI comes from the Haver USECON database (PCU) and the
% euro area CPI comes from the Fagan, Henry and Maestre data base (HICP). 

% XTR
Exports_EA_to_World = log(data(q:tt,5)./data(q:tt,7));
[Exports_EA_to_World_Trend,Exports_EA_to_World_Cycle] = one_sided_hp_filter(Exports_EA_to_World);

% Gross Domestic Product, Billions of Dollars, Quarterly, Seasonally Adjusted Annual Rate
% Transformation:
% Log and HP-Filter to make stationary
% Divided by price level to make real
GDP_US = log(data_US(q:tt,1)./data_US(q:tt,2));
[GDP_US_Trend,GDP_US_Cycle] = one_sided_hp_filter(GDP_US);


% GDP Implicit Price Deflator in United States, Index 2015=100, Quarterly, Seasonally Adjusted
% Transformation:
% Log and First Difference
GDP_Def_US = log(data_US(q-1:tt,2));
GDP_Def_US = (diff(GDP_Def_US));
GDP_Def_US = GDP_Def_US - mean(GDP_Def_US);


% Federal Funds Effective Rate, Percent, Quarterly, Not Seasonally Adjusted
% Transformation:
% Log and Make it Quarterly
FFR = data_US(q:tt,3);
FFR = log(1+FFR/400) - mean(log(1+FFR/400));

% Gross Fixed Capital Formation in United States, US Dollar, Quarterly, Seasonally Adjusted
% Transformation:
% Log and HP-Filter to make stationary
Inve_US = log(data_US(q:tt,4));
[Inve_US_Trend,Inve_US_Cycle] = one_sided_hp_filter(Inve_US);


% % Wage per Head. Calculated as the ratio of compensation of employees, and total employment (WRN = WIN / LNN).
% Wage_US = log((data_US(q:tt,5)./data_US(q:tt,6))./data_US(q:tt,2));
% [Wage_US_Trend,Wage_US_Cycle] = one_sided_hp_filter(Wage_US);
% Wage per Head. Calculated as the ratio of compensation of employees, and total employment (WRN = WIN / LNN).
Wage_US = log((data_US(q:tt,10))./data_US(q:tt,2));
[Wage_US_Trend,Wage_US_Cycle] = one_sided_hp_filter(Wage_US);

% XTIMVA01EZM667S: International Trade: Imports: Value (Goods): Total for the Euro Area (19 Countries), US Dollars, monthly level, Quarterly, Seasonally Adjusted
Imports_US = log(data_US(q:tt,7)./data_US(q:tt,2));
[Imports_US_Trend,Imports_US_Cycle] = one_sided_hp_filter(Imports_US);

% Individual Consumption Expenditure, Millions of euros, Chain linked volume, Calendar and seasonally adjusted data, Reference year 1995.
% Transformation:
% Log and HP-Filter to make stationary
% Divided by price level to make real
Cons_US = log(data_US(q:tt,8)./data_US(q:tt,2));
[Cons_US_Trend,Cons_US_Cycle] = one_sided_hp_filter(Cons_US);

% Zobs = [Euri, GDP_Def_US, GDP_Def, GDP_US_Cycle, GDP_Cycle,  Inve_US_Cycle, Inve_Cycle, Wage_US_Cycle, Wage_Cycle, Exchange_Cycle, FFR ];
% Zobs = [Euri, GDP_Def_US, GDP_Def, GDP_US_Cycle, GDP_Cycle, Inve_US_Cycle, Inve_Cycle, FFR ]; % , Wage_US_Cycle, Wage_Cycle
% Zobs = [Euri, GDP_Def_US, GDP_Def, GDP_US_Cycle, GDP_Cycle, Inve_US_Cycle, Inve_Cycle, FFR, Exchange_Gr,  Wage_US_Cycle, Wage_Cycle ];%,
% Zobs = [Euri, GDP_Def_US, GDP_Def, GDP_US_Cycle, GDP_Cycle,  Inve_US_Cycle, Inve_Cycle, Exchange_Gr, FFR ];
% Zobs = [Euri, GDP_Def_US, GDP_Def, GDP_US_Cycle, GDP_Cycle, FFR ];
% Zobs = [Euri, GDP_Def_US, GDP_Def, FFR ]';


% Zobs = [Euri, GDP_Def_US, GDP_Def, GDP_US_Cycle, GDP_Cycle, Inve_US_Cycle, Inve_Cycle, FFR, Wage_Cycle,  Wage_US_Cycle, Exchange_Gr ];% works well with hessmat290523_Euri_GDP_Def_US_GDP_Def_GDP_US_Cycle_GDP_Cycle_Inve_US_Cycle_Inve_Cycle_FFR_Wage_Cycle_Wage_US_Cycle_Exch_Gr
% Zobs = [Euri, GDP_Def_US, GDP_Def, GDP_US_Cycle, GDP_Cycle, Inve_US_Cycle, Inve_Cycle, FFR, Wage_Cycle,  Wage_US_Cycle, Exchange_n ];
% end













% 
% l=length(FFR);
% 
% DATES = zeros(l,2);
% c = l/4;
% v = [1:c];
% u = repelem(v,4);
% 
% for i = 1:4:l-1;%104
%      for j = 0:3
%         DATES(i+j,1) = (1996+u(i)-1);
%         DATES(i+j,2) = j+1;
%      end
% end
% DATES(105,1) = 2022; DATES(105,2) = 1;
% DATES(106,1) = 2022; DATES(106,2) = 2;


l=length(FFR);

DATES = zeros(l,2);
c = l/4;
v = [1:c];
u = repelem(v,4);

for i = 1:4:l-3;%73
     for j = 0:3
        DATES(i+j,1) = (1999+u(i)-1);
        DATES(i+j,2) = j+1;
     end
end
DATES(73,1) = 2017; DATES(73,2) = 1;
DATES(74,1) = 2017; DATES(74,2) = 2;
DATES(75,1) = 2017; DATES(75,2) = 3;

% this should be used for the estimation in Dynare
% save('Dataset160322_VN','FedFunds', 'RGDP', 'In', 'RWAGE', 'R_V', 'In_V', 'RCons_V', 'Ex_V', 'Inv_V', 'RGDP_V', 'DATES')
% save('Dataset010224','Euri', 'GDP_Cycle', 'GDP_Def', 'Wage_Cycle', 'FFR', 'GDP_US_Cycle', 'GDP_Def_US', 'Wage_US_Cycle', 'Exchange_n', 'DATES')


save('Dataset140524','Euri', 'GDP_Cycle', 'GDP_Def', 'Wage_Cycle', 'FFR', 'GDP_US_Cycle', 'GDP_Def_US', 'Wage_US_Cycle', 'Exchange_r','Exchange_r_diff', 'DATES')




%%

