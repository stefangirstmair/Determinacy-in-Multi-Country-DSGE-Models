function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = PCP_D.static_g1_tt(T, y, x, params);
end
g1 = zeros(105, 105);
g1(1,39)=exp(y(39));
g1(1,42)=(-exp(y(42)));
g1(1,43)=(-exp(y(43)));
g1(1,44)=(-exp(y(44)));
g1(2,41)=exp(y(41));
g1(2,45)=(-exp(y(45)));
g1(2,46)=(-exp(y(46)));
g1(2,47)=(-exp(y(47)));
g1(3,40)=exp(y(40));
g1(3,48)=(-exp(y(48)));
g1(3,49)=(-exp(y(49)));
g1(3,50)=(-exp(y(50)));
g1(4,1)=(-(T(3)*params(8)*(-params(21))*getPowerDeriv(1-params(21)*y(1),T(1),1)));
g1(4,13)=(-(T(2)*exp(y(13))));
g1(4,36)=(-(T(2)*exp(y(36))));
g1(4,42)=exp(y(42));
g1(5,43)=exp(y(43));
g1(5,64)=(-(params(33)*exp(y(64))));
g1(5,65)=(-(params(34)*exp(y(65))));
g1(6,15)=(-(T(5)*exp(y(15))));
g1(6,37)=(-(T(5)*exp(y(37))));
g1(6,54)=T(35);
g1(6,64)=exp(y(64));
g1(6,66)=T(35);
g1(7,15)=(-(exp(y(15))*T(7)));
g1(7,37)=(-(exp(y(37))*T(7)));
g1(7,65)=exp(y(65));
g1(7,67)=(-(T(6)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*y(67),T(1),1)));
g1(8,16)=(-(params(1)*params(61)));
g1(8,33)=(-T(8));
g1(8,54)=(-(T(8)*(-params(22))));
g1(8,68)=1-params(1)*params(61);
g1(9,18)=(-(params(1)*params(61)));
g1(9,33)=(-T(8));
g1(9,54)=(-T(8));
g1(9,69)=1-params(1)*params(61);
g1(10,16)=1-(1-params(61));
g1(10,66)=1-params(61);
g1(10,68)=(-(1-params(61)));
g1(11,18)=1-(1-params(61));
g1(11,67)=1-params(61);
g1(11,69)=(-(1-params(61)));
g1(12,44)=exp(y(44));
g1(12,70)=(-(params(35)*exp(y(70))));
g1(12,71)=(-(params(36)*exp(y(71))));
g1(13,14)=(-(T(9)*exp(y(14))));
g1(13,38)=(-(T(9)*exp(y(38))));
g1(13,55)=T(40);
g1(13,70)=exp(y(70));
g1(13,72)=T(40);
g1(14,14)=(-(exp(y(14))*T(11)));
g1(14,38)=(-(exp(y(38))*T(11)));
g1(14,71)=exp(y(71));
g1(14,73)=(-(T(10)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*y(73),T(1),1)));
g1(15,16)=(-(params(1)*params(61)));
g1(15,33)=(-T(8));
g1(15,55)=(-(T(8)*(-params(22))));
g1(15,74)=1-params(1)*params(61);
g1(16,17)=(-(params(1)*params(61)));
g1(16,33)=(-T(8));
g1(16,55)=(-T(8));
g1(16,75)=1-params(1)*params(61);
g1(17,16)=1-(1-params(61));
g1(17,72)=1-params(61);
g1(17,74)=(-(1-params(61)));
g1(18,17)=1-(1-params(61));
g1(18,73)=1-params(61);
g1(18,75)=(-(1-params(61)));
g1(19,2)=(-(T(6)*params(8)*(-params(21))*getPowerDeriv(1-params(21)*y(2),T(1),1)));
g1(19,15)=(-(exp(y(15))*T(12)));
g1(19,37)=(-(exp(y(37))*T(12)));
g1(19,45)=exp(y(45));
g1(20,46)=exp(y(46));
g1(20,76)=(-(params(37)*exp(y(76))));
g1(20,77)=(-(params(38)*exp(y(77))));
g1(21,13)=(-(exp(y(13))*T(13)));
g1(21,36)=(-(exp(y(36))*T(13)));
g1(21,76)=exp(y(76));
g1(21,78)=(-(T(3)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*y(78),T(1),1)));
g1(22,13)=(-(exp(y(13))*T(14)));
g1(22,36)=(-(exp(y(36))*T(14)));
g1(22,54)=(-(T(3)*T(4)*params(21)*T(36)));
g1(22,77)=exp(y(77));
g1(22,79)=(-(T(3)*T(4)*(-params(21))*T(36)));
g1(23,16)=(-(params(1)*params(64)));
g1(23,34)=(-T(15));
g1(23,54)=T(15);
g1(23,80)=1-params(1)*params(64);
g1(24,18)=(-(params(1)*params(64)));
g1(24,34)=(-T(15));
g1(24,54)=(-(params(22)*T(15)));
g1(24,81)=1-params(1)*params(64);
g1(25,16)=1-(1-params(64));
g1(25,78)=1-params(64);
g1(25,80)=(-(1-params(64)));
g1(26,18)=1-(1-params(64));
g1(26,79)=1-params(64);
g1(26,81)=(-(1-params(64)));
g1(27,47)=exp(y(47));
g1(27,82)=(-(params(39)*exp(y(82))));
g1(27,83)=(-(params(40)*exp(y(83))));
g1(27,84)=(-(params(41)*exp(y(84))));
g1(28,14)=(-(exp(y(14))*T(16)));
g1(28,38)=(-(exp(y(38))*T(16)));
g1(28,82)=exp(y(82));
g1(28,85)=(-(T(10)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*y(85),T(1),1)));
g1(29,14)=(-(exp(y(14))*T(17)));
g1(29,38)=(-(exp(y(38))*T(17)));
g1(29,54)=(-(T(10)*T(4)*params(21)*T(37)));
g1(29,55)=(-(T(10)*T(4)*(-params(21))*T(37)));
g1(29,83)=exp(y(83));
g1(29,86)=(-(T(10)*T(4)*(-params(21))*T(37)));
g1(30,14)=(-(exp(y(14))*T(18)));
g1(30,38)=(-(exp(y(38))*T(18)));
g1(30,55)=T(41);
g1(30,84)=exp(y(84));
g1(30,87)=T(41);
g1(31,17)=(-(params(1)*params(64)));
g1(31,34)=(-T(15));
g1(31,54)=T(15);
g1(31,55)=(-T(15));
g1(31,88)=1-params(1)*params(64);
g1(32,18)=(-(params(1)*params(64)));
g1(32,34)=(-T(15));
g1(32,55)=(-(T(15)*(-params(22))));
g1(32,89)=1-params(1)*params(64);
g1(33,16)=(-(params(1)*params(64)));
g1(33,34)=(-T(15));
g1(33,54)=T(15);
g1(33,55)=(-(T(15)*(-params(22))));
g1(33,90)=1-params(1)*params(64);
g1(34,17)=1-(1-params(64));
g1(34,85)=1-params(64);
g1(34,88)=(-(1-params(64)));
g1(35,18)=1-(1-params(64));
g1(35,86)=1-params(64);
g1(35,89)=(-(1-params(64)));
g1(36,16)=1-(1-params(64));
g1(36,87)=1-params(64);
g1(36,90)=(-(1-params(64)));
g1(37,3)=(-(T(10)*params(8)*(-params(21))*getPowerDeriv(1-params(21)*y(3),T(1),1)));
g1(37,14)=(-(exp(y(14))*T(19)));
g1(37,38)=(-(exp(y(38))*T(19)));
g1(37,48)=exp(y(48));
g1(38,49)=exp(y(49));
g1(38,91)=(-(params(42)*exp(y(91))));
g1(38,92)=(-(params(43)*exp(y(92))));
g1(39,13)=(-(exp(y(13))*T(20)));
g1(39,36)=(-(exp(y(36))*T(20)));
g1(39,91)=exp(y(91));
g1(39,93)=(-(T(3)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*y(93),T(1),1)));
g1(40,13)=(-(exp(y(13))*T(21)));
g1(40,36)=(-(exp(y(36))*T(21)));
g1(40,55)=(-(T(3)*T(4)*params(21)*T(42)));
g1(40,92)=exp(y(92));
g1(40,94)=(-(T(3)*T(4)*(-params(21))*T(42)));
g1(41,16)=(-(params(1)*params(67)));
g1(41,35)=(-T(22));
g1(41,55)=T(22);
g1(41,95)=1-params(1)*params(67);
g1(42,17)=(-(params(1)*params(67)));
g1(42,35)=(-T(22));
g1(42,55)=(-(params(22)*T(22)));
g1(42,96)=1-params(1)*params(67);
g1(43,16)=1-(1-params(67));
g1(43,93)=1-params(67);
g1(43,95)=(-(1-params(67)));
g1(44,17)=1-(1-params(67));
g1(44,94)=1-params(67);
g1(44,96)=(-(1-params(67)));
g1(45,50)=exp(y(50));
g1(45,97)=(-(params(44)*exp(y(97))));
g1(45,98)=(-(params(45)*exp(y(98))));
g1(45,99)=(-(params(46)*exp(y(99))));
g1(46,15)=(-(exp(y(15))*T(23)));
g1(46,37)=(-(exp(y(37))*T(23)));
g1(46,54)=(-(T(6)*T(4)*(-params(21))*T(38)));
g1(46,55)=(-(T(6)*T(4)*params(21)*T(38)));
g1(46,97)=exp(y(97));
g1(46,100)=(-(T(6)*T(4)*(-params(21))*T(38)));
g1(47,15)=(-(exp(y(15))*T(24)));
g1(47,37)=(-(exp(y(37))*T(24)));
g1(47,98)=exp(y(98));
g1(47,101)=(-(T(6)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*y(101),T(1),1)));
g1(48,15)=(-(exp(y(15))*T(25)));
g1(48,37)=(-(exp(y(37))*T(25)));
g1(48,54)=T(39);
g1(48,99)=exp(y(99));
g1(48,102)=T(39);
g1(49,17)=(-(params(1)*params(67)));
g1(49,35)=(-T(22));
g1(49,54)=(-(T(22)*(-params(22))));
g1(49,103)=1-params(1)*params(67);
g1(50,18)=(-(params(1)*params(67)));
g1(50,35)=(-T(22));
g1(50,54)=(-T(22));
g1(50,55)=T(22);
g1(50,104)=1-params(1)*params(67);
g1(51,16)=(-(params(1)*params(67)));
g1(51,35)=(-T(22));
g1(51,54)=(-(T(22)*(-params(22))));
g1(51,55)=T(22);
g1(51,105)=1-params(1)*params(67);
g1(52,17)=1-(1-params(67));
g1(52,100)=1-params(67);
g1(52,103)=(-(1-params(67)));
g1(53,18)=1-(1-params(67));
g1(53,101)=1-params(67);
g1(53,104)=(-(1-params(67)));
g1(54,16)=1-(1-params(67));
g1(54,102)=1-params(67);
g1(54,105)=(-(1-params(67)));
g1(55,1)=params(8);
g1(55,54)=T(4)*(-params(38));
g1(55,55)=T(4)*(-params(43));
g1(55,78)=T(4)*params(37);
g1(55,79)=T(4)*params(38);
g1(55,93)=T(4)*params(42);
g1(55,94)=T(4)*params(43);
g1(56,2)=params(8);
g1(56,54)=T(4)*(params(44)+params(33)+params(46));
g1(56,55)=T(4)*(-params(44));
g1(56,66)=params(33)*T(4);
g1(56,67)=params(34)*T(4);
g1(56,100)=T(4)*params(44);
g1(56,101)=T(4)*params(45);
g1(56,102)=T(4)*params(46);
g1(57,3)=params(8);
g1(57,54)=T(4)*(-params(40));
g1(57,55)=T(4)*(params(41)+params(35)+params(40));
g1(57,72)=T(4)*params(35);
g1(57,73)=T(4)*params(36);
g1(57,85)=T(4)*params(39);
g1(57,86)=T(4)*params(40);
g1(57,87)=T(4)*params(41);
g1(58,10)=1-params(1)*params(61);
g1(58,16)=(-(params(1)*params(61)));
g1(58,33)=(-T(8));
g1(59,11)=1-params(1)*params(64);
g1(59,18)=(-(params(1)*params(64)));
g1(59,34)=(-T(15));
g1(60,12)=1-params(1)*params(67);
g1(60,17)=(-(params(1)*params(67)));
g1(60,35)=(-T(22));
g1(61,1)=1-params(61);
g1(61,10)=(-(1-params(61)));
g1(61,16)=1-(1-params(61));
g1(62,2)=1-params(64);
g1(62,11)=(-(1-params(64)));
g1(62,18)=1-(1-params(64));
g1(63,3)=1-params(67);
g1(63,12)=(-(1-params(67)));
g1(63,17)=1-(1-params(67));
g1(64,5)=1;
g1(64,72)=(-params(35));
g1(64,73)=(-params(36));
g1(65,4)=1;
g1(65,66)=(-params(33));
g1(65,67)=(-params(34));
g1(66,6)=1;
g1(66,78)=(-params(37));
g1(66,79)=(-params(38));
g1(67,7)=1;
g1(67,85)=(-params(39));
g1(67,86)=(-params(40));
g1(67,87)=(-params(41));
g1(68,9)=1;
g1(68,93)=(-params(42));
g1(68,94)=(-params(43));
g1(69,8)=1;
g1(69,100)=(-params(44));
g1(69,101)=(-params(45));
g1(69,102)=(-params(46));
g1(70,13)=T(32)-params(1)*(1+y(19))*T(32)/(1+y(16));
g1(70,16)=(-((-(T(26)*params(1)*(1+y(19))))/((1+y(16))*(1+y(16)))));
g1(70,19)=(-(params(1)*T(26)/(1+y(16))));
g1(71,15)=T(34)-params(1)*(1+y(21))*T(34)/(1+y(18));
g1(71,18)=(-((-(T(27)*params(1)*(1+y(21))))/((1+y(18))*(1+y(18)))));
g1(71,21)=(-(params(1)*T(27)/(1+y(18))));
g1(72,14)=T(33)-params(1)*(1+y(20))*T(33)/(1+y(17));
g1(72,17)=(-((-(T(28)*params(1)*(1+y(20))))/((1+y(17))*(1+y(17)))));
g1(72,20)=(-(params(1)*T(28)/(1+y(17))));
g1(73,15)=exp(y(15));
g1(73,16)=(-(exp(y(54))*(1+y(22))*y(25)))/((1+y(16))*(1+y(16)));
g1(73,22)=exp(y(54))*y(25)/(1+y(16));
g1(73,25)=exp(y(54))*(1+y(22))/(1+y(16))-exp(y(54));
g1(73,28)=(-exp(y(28)+y(31)));
g1(73,31)=(-exp(y(28)+y(31)));
g1(73,52)=(-1);
g1(73,54)=exp(y(54))*(1+y(22))*y(25)/(1+y(16))-exp(y(54))*y(25);
g1(74,14)=exp(y(14));
g1(74,16)=(-(exp(y(55))*(1+y(23))*y(26)))/((1+y(16))*(1+y(16)));
g1(74,23)=exp(y(55))*y(26)/(1+y(16));
g1(74,26)=exp(y(55))*(1+y(23))/(1+y(16))-exp(y(55));
g1(74,29)=(-exp(y(29)+y(32)));
g1(74,32)=(-exp(y(29)+y(32)));
g1(74,53)=(-1);
g1(74,55)=exp(y(55))*(1+y(23))*y(26)/(1+y(16))-exp(y(55))*y(26);
g1(75,13)=(-(params(5)*T(29)));
g1(75,16)=1-params(1);
g1(75,27)=T(29);
g1(75,30)=(-(params(6)*T(29)));
g1(76,15)=(-(params(5)*T(30)));
g1(76,18)=1-params(1);
g1(76,28)=T(30);
g1(76,31)=(-(params(6)*T(30)));
g1(77,14)=(-(params(5)*T(31)));
g1(77,17)=1-params(1);
g1(77,29)=T(31);
g1(77,32)=(-(params(6)*T(31)));
g1(78,19)=(-(params(20)/(1+y(19)*params(20))));
g1(78,27)=(-(1-params(2)));
g1(78,33)=1;
g1(78,59)=1;
g1(79,22)=(-(params(20)/(1+y(22)*params(20))));
g1(79,28)=(-(1-params(2)));
g1(79,34)=1;
g1(79,60)=1;
g1(80,23)=(-(params(20)/(1+y(23)*params(20))));
g1(80,29)=(-(1-params(2)));
g1(80,35)=1;
g1(80,61)=1;
g1(81,27)=(-(exp(y(27))/exp(y(33))));
g1(81,30)=(-(exp(y(39))*(1-params(2))*exp(y(30))))/(exp(y(30))*exp(y(30)));
g1(81,33)=(-((-(exp(y(27))*exp(y(33))))/(exp(y(33))*exp(y(33)))));
g1(81,39)=exp(y(39))*(1-params(2))/exp(y(30));
g1(82,28)=(-(exp(y(28))/exp(y(34))));
g1(82,31)=(-(exp(y(41))*(1-params(2))*exp(y(31))))/(exp(y(31))*exp(y(31)));
g1(82,34)=(-((-(exp(y(28))*exp(y(34))))/(exp(y(34))*exp(y(34)))));
g1(82,41)=exp(y(41))*(1-params(2))/exp(y(31));
g1(83,29)=(-(exp(y(29))/exp(y(35))));
g1(83,32)=(-(exp(y(40))*(1-params(2))*exp(y(32))))/(exp(y(32))*exp(y(32)));
g1(83,35)=(-((-(exp(y(29))*exp(y(35))))/(exp(y(35))*exp(y(35)))));
g1(83,40)=exp(y(40))*(1-params(2))/exp(y(32));
g1(84,33)=(-((-exp(y(33)))/(exp(y(33))*exp(y(33)))));
g1(84,36)=(-(exp(y(36))*exp(y(39))*params(2)))/(exp(y(36))*exp(y(36)));
g1(84,39)=exp(y(39))*params(2)/exp(y(36));
g1(85,34)=(-((-exp(y(34)))/(exp(y(34))*exp(y(34)))));
g1(85,37)=(-(exp(y(37))*exp(y(41))*params(2)))/(exp(y(37))*exp(y(37)));
g1(85,41)=exp(y(41))*params(2)/exp(y(37));
g1(86,35)=(-((-exp(y(35)))/(exp(y(35))*exp(y(35)))));
g1(86,38)=(-(exp(y(38))*exp(y(40))*params(2)))/(exp(y(38))*exp(y(38)));
g1(86,40)=exp(y(40))*params(2)/exp(y(38));
g1(87,1)=(-exp(y(42)+y(1)));
g1(87,4)=(-exp(y(43)+y(4)));
g1(87,5)=(-exp(y(44)+y(5)));
g1(87,33)=exp(y(39)+y(33));
g1(87,39)=exp(y(39)+y(33));
g1(87,42)=(-exp(y(42)+y(1)));
g1(87,43)=(-exp(y(43)+y(4)));
g1(87,44)=(-exp(y(44)+y(5)));
g1(87,51)=1;
g1(88,2)=(-exp(y(45)+y(2)));
g1(88,6)=(-exp(y(46)+y(54)+y(6)));
g1(88,7)=(-exp(y(47)+y(54)+y(7)));
g1(88,34)=exp(y(41)+y(34));
g1(88,41)=exp(y(41)+y(34));
g1(88,45)=(-exp(y(45)+y(2)));
g1(88,46)=(-exp(y(46)+y(54)+y(6)));
g1(88,47)=(-exp(y(47)+y(54)+y(7)));
g1(88,52)=1;
g1(88,54)=(-(exp(y(46)+y(54)+y(6))+exp(y(47)+y(54)+y(7))));
g1(89,3)=(-exp(y(48)+y(3)));
g1(89,8)=(-exp(y(50)+y(55)+y(8)));
g1(89,9)=(-exp(y(49)+y(55)+y(9)));
g1(89,35)=exp(y(40)+y(35));
g1(89,40)=exp(y(40)+y(35));
g1(89,48)=(-exp(y(48)+y(3)));
g1(89,49)=(-exp(y(49)+y(55)+y(9)));
g1(89,50)=(-exp(y(50)+y(55)+y(8)));
g1(89,53)=1;
g1(89,55)=(-(exp(y(49)+y(55)+y(9))+exp(y(50)+y(55)+y(8))));
g1(90,19)=(-1);
g1(90,22)=1;
g1(90,25)=(-(params(10)*exp(y(25)-params(11))));
g1(91,19)=(-1);
g1(91,23)=1;
g1(91,26)=(-(params(10)*exp(y(26)-params(11))));
g1(92,16)=(-((1-params(17))*params(60)));
g1(92,19)=1-params(17);
g1(92,39)=(-((1-params(17))*params(30)));
g1(92,56)=1;
g1(93,18)=(-((1-params(17))*params(60)));
g1(93,21)=1-params(17);
g1(93,41)=(-((1-params(17))*params(30)));
g1(93,57)=1;
g1(94,17)=(-((1-params(17))*params(60)));
g1(94,20)=1-params(17);
g1(94,40)=(-((1-params(17))*params(30)));
g1(94,58)=1;
g1(95,16)=1;
g1(95,18)=(-1);
g1(95,21)=1;
g1(95,22)=(-1);
g1(95,62)=(-1);
g1(96,16)=1;
g1(96,17)=(-1);
g1(96,20)=1;
g1(96,23)=(-1);
g1(96,63)=(-1);
g1(97,24)=1;
g1(97,25)=1;
g1(97,26)=1;
g1(98,56)=1-params(18);
g1(99,57)=1-params(18);
g1(100,58)=1-params(18);
g1(101,62)=1-params(28);
g1(102,63)=1-params(28);
g1(103,59)=1-params(26);
g1(104,60)=1-params(26);
g1(105,61)=1-params(26);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
