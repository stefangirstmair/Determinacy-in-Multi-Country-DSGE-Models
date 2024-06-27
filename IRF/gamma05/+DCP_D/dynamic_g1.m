function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = DCP_D.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(105, 181);
g1(1,77)=exp(y(77));
g1(1,80)=(-exp(y(80)));
g1(1,81)=(-exp(y(81)));
g1(1,82)=(-exp(y(82)));
g1(2,79)=exp(y(79));
g1(2,83)=(-exp(y(83)));
g1(2,84)=(-exp(y(84)));
g1(2,85)=(-exp(y(85)));
g1(3,78)=exp(y(78));
g1(3,86)=(-exp(y(86)));
g1(3,87)=(-exp(y(87)));
g1(3,88)=(-exp(y(88)));
g1(4,39)=(-(T(3)*params(8)*(-params(21))*getPowerDeriv(1-params(21)*y(39),T(1),1)));
g1(4,51)=(-(T(2)*exp(y(51))));
g1(4,74)=(-(T(2)*exp(y(74))));
g1(4,80)=exp(y(80));
g1(5,81)=exp(y(81));
g1(5,102)=(-(params(33)*exp(y(102))));
g1(5,103)=(-(params(34)*exp(y(103))));
g1(6,53)=(-(T(5)*exp(y(53))));
g1(6,75)=(-(T(5)*exp(y(75))));
g1(6,92)=T(39);
g1(6,102)=exp(y(102));
g1(6,104)=T(39);
g1(7,53)=(-(exp(y(53))*T(7)));
g1(7,75)=(-(exp(y(75))*T(7)));
g1(7,103)=exp(y(103));
g1(7,105)=(-(T(6)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*y(105),T(1),1)));
g1(8,150)=T(32);
g1(8,71)=(-T(8));
g1(8,92)=(-(T(8)*(-params(22))));
g1(8,106)=1;
g1(8,158)=T(32);
g1(9,152)=T(32);
g1(9,71)=(-T(8));
g1(9,92)=(-T(8));
g1(9,107)=1;
g1(9,159)=T(32);
g1(10,54)=1-(1-params(14));
g1(10,25)=1-params(14)-1;
g1(10,104)=1;
g1(10,106)=(-(1-params(14)));
g1(11,56)=1-(1-params(14));
g1(11,26)=1-params(14)-1;
g1(11,105)=1;
g1(11,107)=(-(1-params(14)));
g1(12,82)=exp(y(82));
g1(12,108)=(-(params(35)*exp(y(108))));
g1(12,109)=(-(params(36)*exp(y(109))));
g1(13,52)=(-(T(9)*exp(y(52))));
g1(13,76)=(-(T(9)*exp(y(76))));
g1(13,93)=T(44);
g1(13,108)=exp(y(108));
g1(13,110)=T(44);
g1(14,52)=(-(exp(y(52))*T(11)));
g1(14,76)=(-(exp(y(76))*T(11)));
g1(14,109)=exp(y(109));
g1(14,111)=(-(T(10)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*y(111),T(1),1)));
g1(15,150)=T(32);
g1(15,71)=(-T(8));
g1(15,93)=(-(T(8)*(-params(22))));
g1(15,112)=1;
g1(15,160)=T(32);
g1(16,151)=T(32);
g1(16,71)=(-T(8));
g1(16,93)=(-T(8));
g1(16,113)=1;
g1(16,161)=T(32);
g1(17,54)=1-(1-params(14));
g1(17,27)=1-params(14)-1;
g1(17,110)=1;
g1(17,112)=(-(1-params(14)));
g1(18,55)=1-(1-params(14));
g1(18,28)=1-params(14)-1;
g1(18,111)=1;
g1(18,113)=(-(1-params(14)));
g1(19,40)=(-(T(6)*params(8)*(-params(21))*getPowerDeriv(1-params(21)*y(40),T(1),1)));
g1(19,53)=(-(exp(y(53))*T(12)));
g1(19,75)=(-(exp(y(75))*T(12)));
g1(19,83)=exp(y(83));
g1(20,84)=exp(y(84));
g1(20,114)=(-(params(37)*exp(y(114))));
g1(20,115)=(-(params(38)*exp(y(115))));
g1(21,51)=(-(exp(y(51))*T(13)));
g1(21,74)=(-(exp(y(74))*T(13)));
g1(21,114)=exp(y(114));
g1(21,116)=(-(T(3)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*y(116),T(1),1)));
g1(22,51)=(-(exp(y(51))*T(14)));
g1(22,74)=(-(exp(y(74))*T(14)));
g1(22,92)=(-(T(3)*T(4)*params(21)*T(40)));
g1(22,115)=exp(y(115));
g1(22,117)=(-(T(3)*T(4)*(-params(21))*T(40)));
g1(23,150)=T(32);
g1(23,72)=(-T(8));
g1(23,92)=T(8);
g1(23,118)=1;
g1(23,162)=T(32);
g1(24,152)=T(32);
g1(24,72)=(-T(8));
g1(24,92)=(-(params(22)*T(8)));
g1(24,119)=1;
g1(24,163)=T(32);
g1(25,54)=1-(1-params(14));
g1(25,29)=1-params(14)-1;
g1(25,116)=1;
g1(25,118)=(-(1-params(14)));
g1(26,56)=1-(1-params(14));
g1(26,30)=1-params(14)-1;
g1(26,117)=1;
g1(26,119)=(-(1-params(14)));
g1(27,85)=exp(y(85));
g1(27,120)=(-(params(39)*exp(y(120))));
g1(27,121)=(-(params(40)*exp(y(121))));
g1(27,122)=(-(params(41)*exp(y(122))));
g1(28,52)=(-(exp(y(52))*T(15)));
g1(28,76)=(-(exp(y(76))*T(15)));
g1(28,120)=exp(y(120));
g1(28,123)=(-(T(10)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*y(123),T(1),1)));
g1(29,52)=(-(exp(y(52))*T(16)));
g1(29,76)=(-(exp(y(76))*T(16)));
g1(29,92)=(-(T(10)*T(4)*params(21)*T(41)));
g1(29,93)=(-(T(10)*T(4)*(-params(21))*T(41)));
g1(29,121)=exp(y(121));
g1(29,124)=(-(T(10)*T(4)*(-params(21))*T(41)));
g1(30,52)=(-(exp(y(52))*T(17)));
g1(30,76)=(-(exp(y(76))*T(17)));
g1(30,93)=T(45);
g1(30,122)=exp(y(122));
g1(30,125)=T(45);
g1(31,151)=T(32);
g1(31,72)=(-T(8));
g1(31,92)=T(8);
g1(31,93)=(-T(8));
g1(31,126)=1;
g1(31,164)=T(32);
g1(32,152)=T(32);
g1(32,72)=(-T(8));
g1(32,93)=(-(T(8)*(-params(22))));
g1(32,127)=1;
g1(32,165)=T(32);
g1(33,150)=T(32);
g1(33,72)=(-T(8));
g1(33,92)=T(8);
g1(33,93)=(-(T(8)*(-params(22))));
g1(33,128)=1;
g1(33,166)=T(32);
g1(34,55)=1-(1-params(14));
g1(34,31)=1-params(14)-1;
g1(34,123)=1;
g1(34,126)=(-(1-params(14)));
g1(35,56)=1-(1-params(14));
g1(35,32)=1-params(14)-1;
g1(35,124)=1;
g1(35,127)=(-(1-params(14)));
g1(36,54)=1-(1-params(14));
g1(36,33)=1-params(14)-1;
g1(36,125)=1;
g1(36,128)=(-(1-params(14)));
g1(37,41)=(-(T(10)*params(8)*(-params(21))*getPowerDeriv(1-params(21)*y(41),T(1),1)));
g1(37,52)=(-(exp(y(52))*T(18)));
g1(37,76)=(-(exp(y(76))*T(18)));
g1(37,86)=exp(y(86));
g1(38,87)=exp(y(87));
g1(38,129)=(-(params(42)*exp(y(129))));
g1(38,130)=(-(params(43)*exp(y(130))));
g1(39,51)=(-(exp(y(51))*T(19)));
g1(39,74)=(-(exp(y(74))*T(19)));
g1(39,129)=exp(y(129));
g1(39,131)=(-(T(3)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*y(131),T(1),1)));
g1(40,51)=(-(exp(y(51))*T(20)));
g1(40,74)=(-(exp(y(74))*T(20)));
g1(40,93)=(-(T(3)*T(4)*params(21)*T(46)));
g1(40,130)=exp(y(130));
g1(40,132)=(-(T(3)*T(4)*(-params(21))*T(46)));
g1(41,150)=T(32);
g1(41,73)=(-T(8));
g1(41,93)=T(8);
g1(41,133)=1;
g1(41,167)=T(32);
g1(42,151)=T(32);
g1(42,73)=(-T(8));
g1(42,93)=(-(params(22)*T(8)));
g1(42,134)=1;
g1(42,168)=T(32);
g1(43,54)=1-(1-params(14));
g1(43,34)=1-params(14)-1;
g1(43,131)=1;
g1(43,133)=(-(1-params(14)));
g1(44,55)=1-(1-params(14));
g1(44,35)=1-params(14)-1;
g1(44,132)=1;
g1(44,134)=(-(1-params(14)));
g1(45,88)=exp(y(88));
g1(45,135)=(-(params(44)*exp(y(135))));
g1(45,136)=(-(params(45)*exp(y(136))));
g1(45,137)=(-(params(46)*exp(y(137))));
g1(46,53)=(-(exp(y(53))*T(21)));
g1(46,75)=(-(exp(y(75))*T(21)));
g1(46,92)=(-(T(6)*T(4)*(-params(21))*T(42)));
g1(46,93)=(-(T(6)*T(4)*params(21)*T(42)));
g1(46,135)=exp(y(135));
g1(46,138)=(-(T(6)*T(4)*(-params(21))*T(42)));
g1(47,53)=(-(exp(y(53))*T(22)));
g1(47,75)=(-(exp(y(75))*T(22)));
g1(47,136)=exp(y(136));
g1(47,139)=(-(T(6)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*y(139),T(1),1)));
g1(48,53)=(-(exp(y(53))*T(23)));
g1(48,75)=(-(exp(y(75))*T(23)));
g1(48,92)=T(43);
g1(48,137)=exp(y(137));
g1(48,140)=T(43);
g1(49,151)=T(32);
g1(49,73)=(-T(8));
g1(49,92)=(-(T(8)*(-params(22))));
g1(49,141)=1;
g1(49,169)=T(32);
g1(50,152)=T(32);
g1(50,73)=(-T(8));
g1(50,92)=(-T(8));
g1(50,93)=T(8);
g1(50,142)=1;
g1(50,170)=T(32);
g1(51,150)=T(32);
g1(51,73)=(-T(8));
g1(51,92)=(-(T(8)*(-params(22))));
g1(51,93)=T(8);
g1(51,143)=1;
g1(51,171)=T(32);
g1(52,55)=1-(1-params(14));
g1(52,36)=1-params(14)-1;
g1(52,138)=1;
g1(52,141)=(-(1-params(14)));
g1(53,56)=1-(1-params(14));
g1(53,37)=1-params(14)-1;
g1(53,139)=1;
g1(53,142)=(-(1-params(14)));
g1(54,54)=1-(1-params(14));
g1(54,38)=1-params(14)-1;
g1(54,140)=1;
g1(54,143)=(-(1-params(14)));
g1(55,39)=params(8);
g1(55,92)=T(4)*(-params(38));
g1(55,93)=T(4)*(-params(43));
g1(55,116)=T(4)*params(37);
g1(55,117)=T(4)*params(38);
g1(55,131)=T(4)*params(42);
g1(55,132)=T(4)*params(43);
g1(56,40)=params(8);
g1(56,92)=T(4)*(params(44)+params(33)+params(46));
g1(56,93)=T(4)*(-params(44));
g1(56,104)=params(33)*T(4);
g1(56,105)=params(34)*T(4);
g1(56,138)=T(4)*params(44);
g1(56,139)=T(4)*params(45);
g1(56,140)=T(4)*params(46);
g1(57,41)=params(8);
g1(57,92)=T(4)*(-params(40));
g1(57,93)=T(4)*(params(41)+params(35)+params(40));
g1(57,110)=T(4)*params(35);
g1(57,111)=T(4)*params(36);
g1(57,123)=T(4)*params(39);
g1(57,124)=T(4)*params(40);
g1(57,125)=T(4)*params(41);
g1(58,48)=1;
g1(58,144)=T(32);
g1(58,150)=T(32);
g1(58,71)=(-T(8));
g1(59,49)=1;
g1(59,145)=T(32);
g1(59,152)=T(32);
g1(59,72)=(-T(8));
g1(60,50)=1;
g1(60,146)=T(32);
g1(60,151)=T(32);
g1(60,73)=(-T(8));
g1(61,1)=1-params(14)-1;
g1(61,39)=1;
g1(61,48)=(-(1-params(14)));
g1(61,54)=1-(1-params(14));
g1(62,2)=1-params(14)-1;
g1(62,40)=1;
g1(62,49)=(-(1-params(14)));
g1(62,56)=1-(1-params(14));
g1(63,3)=1-params(14)-1;
g1(63,41)=1;
g1(63,50)=(-(1-params(14)));
g1(63,55)=1-(1-params(14));
g1(64,43)=1;
g1(64,110)=(-params(35));
g1(64,111)=(-params(36));
g1(65,47)=1;
g1(65,131)=(-params(42));
g1(65,132)=(-params(43));
g1(66,42)=1;
g1(66,104)=(-params(33));
g1(66,105)=(-params(34));
g1(67,44)=1;
g1(67,116)=(-params(37));
g1(67,117)=(-params(38));
g1(68,45)=1;
g1(68,123)=(-params(39));
g1(68,124)=(-params(40));
g1(68,125)=(-params(41));
g1(69,46)=1;
g1(69,138)=(-params(44));
g1(69,139)=(-params(45));
g1(69,140)=(-params(46));
g1(70,39)=(-exp(y(80)+y(39)));
g1(70,71)=exp(y(77)+y(71));
g1(70,77)=exp(y(77)+y(71));
g1(70,80)=(-exp(y(80)+y(39)));
g1(70,89)=1;
g1(70,93)=(-(params(36)*exp(y(109)+y(93)+y(111))));
g1(70,102)=(-(params(33)*exp(y(102)+y(104))));
g1(70,103)=(-(params(34)*exp(y(103)+y(105))));
g1(70,104)=(-(params(33)*exp(y(102)+y(104))));
g1(70,105)=(-(params(34)*exp(y(103)+y(105))));
g1(70,108)=(-(params(35)*exp(y(108)+y(110))));
g1(70,109)=(-(params(36)*exp(y(109)+y(93)+y(111))));
g1(70,110)=(-(params(35)*exp(y(108)+y(110))));
g1(70,111)=(-(params(36)*exp(y(109)+y(93)+y(111))));
g1(71,40)=(-exp(y(83)+y(40)));
g1(71,72)=exp(y(79)+y(72));
g1(71,79)=exp(y(79)+y(72));
g1(71,83)=(-exp(y(83)+y(40)));
g1(71,90)=1;
g1(71,92)=(-(params(39)*exp(y(120)+y(92)+y(123)-y(93))+params(37)*exp(y(114)+y(92)+y(116))+params(41)*exp(y(122)+y(92)+y(125))));
g1(71,93)=(-(params(39)*(-exp(y(120)+y(92)+y(123)-y(93)))));
g1(71,114)=(-(params(37)*exp(y(114)+y(92)+y(116))));
g1(71,115)=(-(params(38)*exp(y(115)+y(117))));
g1(71,116)=(-(params(37)*exp(y(114)+y(92)+y(116))));
g1(71,117)=(-(params(38)*exp(y(115)+y(117))));
g1(71,120)=(-(params(39)*exp(y(120)+y(92)+y(123)-y(93))));
g1(71,121)=(-(params(40)*exp(y(121)+y(124))));
g1(71,122)=(-(params(41)*exp(y(122)+y(92)+y(125))));
g1(71,123)=(-(params(39)*exp(y(120)+y(92)+y(123)-y(93))));
g1(71,124)=(-(params(40)*exp(y(121)+y(124))));
g1(71,125)=(-(params(41)*exp(y(122)+y(92)+y(125))));
g1(72,41)=(-exp(y(86)+y(41)));
g1(72,73)=exp(y(78)+y(73));
g1(72,78)=exp(y(78)+y(73));
g1(72,86)=(-exp(y(86)+y(41)));
g1(72,91)=1;
g1(72,92)=(-(params(45)*(-exp(y(136)+y(93)+y(139)-y(92)))));
g1(72,93)=(-(params(45)*exp(y(136)+y(93)+y(139)-y(92))+params(42)*exp(y(129)+y(93)+y(131))+params(46)*exp(y(137)+y(93)+y(140))));
g1(72,129)=(-(params(42)*exp(y(129)+y(93)+y(131))));
g1(72,130)=(-(params(43)*exp(y(130)+y(132))));
g1(72,131)=(-(params(42)*exp(y(129)+y(93)+y(131))));
g1(72,132)=(-(params(43)*exp(y(130)+y(132))));
g1(72,135)=(-(params(44)*exp(y(135)+y(138))));
g1(72,136)=(-(params(45)*exp(y(136)+y(93)+y(139)-y(92))));
g1(72,137)=(-(params(46)*exp(y(137)+y(93)+y(140))));
g1(72,138)=(-(params(44)*exp(y(135)+y(138))));
g1(72,139)=(-(params(45)*exp(y(136)+y(93)+y(139)-y(92))));
g1(72,140)=(-(params(46)*exp(y(137)+y(93)+y(140))));
g1(73,4)=exp(y(51)-params(32)*y(4))*(-params(32))*T(33);
g1(73,51)=exp(y(51)-params(32)*y(4))*T(33)-params(1)*(1+y(57))*exp(y(147)-y(51)*params(32))*(-params(32))*T(34)/(1+y(150));
g1(73,147)=(-(params(1)*(1+y(57))*exp(y(147)-y(51)*params(32))*T(34)/(1+y(150))));
g1(73,150)=(-((-T(25))/((1+y(150))*(1+y(150)))));
g1(73,57)=(-(params(1)*T(24)/(1+y(150))));
g1(74,6)=exp(y(53)-params(32)*y(6))*(-params(32))*T(37);
g1(74,53)=exp(y(53)-params(32)*y(6))*T(37)-params(1)*(1+y(59))*exp(y(149)-y(53)*params(32))*(-params(32))*T(38)/(1+y(152));
g1(74,149)=(-(params(1)*(1+y(59))*exp(y(149)-y(53)*params(32))*T(38)/(1+y(152))));
g1(74,152)=(-((-T(27))/((1+y(152))*(1+y(152)))));
g1(74,59)=(-(params(1)*T(26)/(1+y(152))));
g1(75,5)=exp(y(52)-params(32)*y(5))*(-params(32))*T(35);
g1(75,52)=exp(y(52)-params(32)*y(5))*T(35)-params(1)*(1+y(58))*exp(y(148)-y(52)*params(32))*(-params(32))*T(36)/(1+y(151));
g1(75,148)=(-(params(1)*(1+y(58))*exp(y(148)-y(52)*params(32))*T(36)/(1+y(151))));
g1(75,151)=(-((-T(29))/((1+y(151))*(1+y(151)))));
g1(75,58)=(-(params(1)*T(28)/(1+y(151))));
g1(76,53)=exp(y(53));
g1(76,54)=(-(exp(y(92))*(1+y(10))*y(12)))/((1+y(54))*(1+y(54)));
g1(76,10)=exp(y(92))*y(12)/(1+y(54));
g1(76,12)=exp(y(92))*(1+y(10))/(1+y(54));
g1(76,63)=(-exp(y(92)));
g1(76,66)=(-exp(y(66)+y(69)));
g1(76,69)=(-exp(y(66)+y(69)));
g1(76,90)=(-1);
g1(76,92)=exp(y(92))*(1+y(10))*y(12)/(1+y(54))-exp(y(92))*y(63);
g1(77,52)=exp(y(52));
g1(77,54)=(-(exp(y(93))*(1+y(11))*y(13)))/((1+y(54))*(1+y(54)));
g1(77,11)=exp(y(93))*y(13)/(1+y(54));
g1(77,13)=exp(y(93))*(1+y(11))/(1+y(54));
g1(77,64)=(-exp(y(93)));
g1(77,67)=(-exp(y(67)+y(70)));
g1(77,70)=(-exp(y(67)+y(70)));
g1(77,91)=(-1);
g1(77,93)=exp(y(93))*(1+y(11))*y(13)/(1+y(54))-exp(y(93))*y(64);
g1(78,4)=(-(T(30)*T(31)*(-params(32))));
g1(78,51)=(-(T(30)*T(31)));
g1(78,54)=1;
g1(78,150)=(-params(1));
g1(78,14)=(-1);
g1(78,65)=1-((-params(1))-T(30));
g1(78,153)=(-params(1));
g1(78,68)=(-(params(6)*T(30)));
g1(79,6)=(-(T(30)*T(31)*(-params(32))));
g1(79,53)=(-(T(30)*T(31)));
g1(79,56)=1;
g1(79,152)=(-params(1));
g1(79,15)=(-1);
g1(79,66)=1-((-params(1))-T(30));
g1(79,154)=(-params(1));
g1(79,69)=(-(params(6)*T(30)));
g1(80,5)=(-(T(30)*T(31)*(-params(32))));
g1(80,52)=(-(T(30)*T(31)));
g1(80,55)=1;
g1(80,151)=(-params(1));
g1(80,16)=(-1);
g1(80,67)=1-((-params(1))-T(30));
g1(80,155)=(-params(1));
g1(80,70)=(-(params(6)*T(30)));
g1(81,7)=(-(params(20)/(1+params(20)*y(7))));
g1(81,65)=(-(1-params(2)));
g1(81,71)=1;
g1(81,97)=1;
g1(82,10)=(-(params(20)/(1+y(10)*params(20))));
g1(82,66)=(-(1-params(2)));
g1(82,72)=1;
g1(82,98)=1;
g1(83,11)=(-(params(20)/(1+y(11)*params(20))));
g1(83,67)=(-(1-params(2)));
g1(83,73)=1;
g1(83,99)=1;
g1(84,65)=(-(exp(y(65))/exp(y(71))));
g1(84,68)=(-(exp(y(77))*(1-params(2))*exp(y(68))))/(exp(y(68))*exp(y(68)));
g1(84,71)=(-((-(exp(y(65))*exp(y(71))))/(exp(y(71))*exp(y(71)))));
g1(84,77)=exp(y(77))*(1-params(2))/exp(y(68));
g1(85,66)=(-(exp(y(66))/exp(y(72))));
g1(85,69)=(-(exp(y(79))*(1-params(2))*exp(y(69))))/(exp(y(69))*exp(y(69)));
g1(85,72)=(-((-(exp(y(66))*exp(y(72))))/(exp(y(72))*exp(y(72)))));
g1(85,79)=exp(y(79))*(1-params(2))/exp(y(69));
g1(86,67)=(-(exp(y(67))/exp(y(73))));
g1(86,70)=(-(exp(y(78))*(1-params(2))*exp(y(70))))/(exp(y(70))*exp(y(70)));
g1(86,73)=(-((-(exp(y(67))*exp(y(73))))/(exp(y(73))*exp(y(73)))));
g1(86,78)=exp(y(78))*(1-params(2))/exp(y(70));
g1(87,71)=(-((-exp(y(71)))/(exp(y(71))*exp(y(71)))));
g1(87,74)=(-(exp(y(74))*exp(y(77))*params(2)))/(exp(y(74))*exp(y(74)));
g1(87,77)=exp(y(77))*params(2)/exp(y(74));
g1(88,72)=(-((-exp(y(72)))/(exp(y(72))*exp(y(72)))));
g1(88,75)=(-(exp(y(75))*exp(y(79))*params(2)))/(exp(y(75))*exp(y(75)));
g1(88,79)=exp(y(79))*params(2)/exp(y(75));
g1(89,73)=(-((-exp(y(73)))/(exp(y(73))*exp(y(73)))));
g1(89,76)=(-(exp(y(76))*exp(y(78))*params(2)))/(exp(y(76))*exp(y(76)));
g1(89,78)=exp(y(78))*params(2)/exp(y(76));
g1(90,57)=(-1);
g1(90,60)=1;
g1(90,63)=(-(params(10)*exp(y(63)-params(11))));
g1(90,176)=(-1);
g1(91,57)=(-1);
g1(91,61)=1;
g1(91,64)=(-(params(10)*exp(y(64)-params(11))));
g1(91,175)=(-1);
g1(92,54)=(-((1-params(53))*params(47)));
g1(92,7)=(-params(53));
g1(92,57)=1;
g1(92,77)=(-((1-params(53))*params(50)));
g1(92,94)=(-1);
g1(93,56)=(-((1-params(54))*params(48)));
g1(93,9)=(-params(54));
g1(93,59)=1;
g1(93,79)=(-((1-params(54))*params(51)));
g1(93,95)=(-1);
g1(94,55)=(-((1-params(17))*params(49)));
g1(94,8)=(-params(17));
g1(94,58)=1;
g1(94,78)=(-((1-params(17))*params(52)));
g1(94,96)=(-1);
g1(95,150)=1;
g1(95,152)=(-1);
g1(95,59)=1;
g1(95,60)=(-1);
g1(95,92)=1;
g1(95,156)=(-1);
g1(95,100)=(-1);
g1(96,150)=1;
g1(96,151)=(-1);
g1(96,58)=1;
g1(96,61)=(-1);
g1(96,93)=1;
g1(96,157)=(-1);
g1(96,101)=(-1);
g1(97,62)=1;
g1(97,63)=1;
g1(97,64)=1;
g1(98,17)=(-params(18));
g1(98,94)=1;
g1(98,172)=(-1);
g1(99,18)=(-params(18));
g1(99,95)=1;
g1(99,173)=(-1);
g1(100,19)=(-params(18));
g1(100,96)=1;
g1(100,174)=(-1);
g1(101,23)=(-params(28));
g1(101,100)=1;
g1(101,180)=(-1);
g1(102,24)=(-params(28));
g1(102,101)=1;
g1(102,181)=(-1);
g1(103,20)=(-params(26));
g1(103,97)=1;
g1(103,177)=(-1);
g1(104,21)=(-params(26));
g1(104,98)=1;
g1(104,178)=(-1);
g1(105,22)=(-params(26));
g1(105,99)=1;
g1(105,179)=(-1);

end
