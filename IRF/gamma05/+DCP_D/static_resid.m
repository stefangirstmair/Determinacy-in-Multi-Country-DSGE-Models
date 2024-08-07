function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
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
%   residual
%

if T_flag
    T = DCP_D.static_resid_tt(T, y, x, params);
end
residual = zeros(105, 1);
lhs = exp(y(39));
rhs = exp(y(42))+exp(y(43))+exp(y(44));
residual(1) = lhs - rhs;
lhs = exp(y(41));
rhs = exp(y(45))+exp(y(46))+exp(y(47));
residual(2) = lhs - rhs;
lhs = exp(y(40));
rhs = exp(y(48))+exp(y(50))+exp(y(49));
residual(3) = lhs - rhs;
lhs = exp(y(42));
rhs = T(2)*T(3);
residual(4) = lhs - rhs;
lhs = exp(y(43));
rhs = params(33)*exp(y(64))+params(34)*exp(y(65));
residual(5) = lhs - rhs;
lhs = exp(y(64));
rhs = T(5)*T(6);
residual(6) = lhs - rhs;
lhs = exp(y(65));
rhs = T(6)*T(7);
residual(7) = lhs - rhs;
lhs = y(68);
rhs = params(1)*params(14)*(y(68)+y(16))+T(8)*(y(33)-y(54)*params(22)+params(23));
residual(8) = lhs - rhs;
lhs = y(69);
rhs = params(1)*params(14)*(y(69)+y(18))+T(8)*(params(23)+y(54)+y(33));
residual(9) = lhs - rhs;
lhs = y(16);
rhs = (1-params(14))*(y(16)+y(68)-y(66));
residual(10) = lhs - rhs;
lhs = y(18);
rhs = (1-params(14))*(y(18)+y(69)-y(67));
residual(11) = lhs - rhs;
lhs = exp(y(44));
rhs = params(35)*exp(y(70))+params(36)*exp(y(71));
residual(12) = lhs - rhs;
lhs = exp(y(70));
rhs = T(9)*T(10);
residual(13) = lhs - rhs;
lhs = exp(y(71));
rhs = T(10)*T(11);
residual(14) = lhs - rhs;
lhs = y(74);
rhs = params(1)*params(14)*(y(16)+y(74))+T(8)*(params(23)+y(33)-params(22)*y(55));
residual(15) = lhs - rhs;
lhs = y(75);
rhs = params(1)*params(14)*(y(75)+y(17))+T(8)*(params(23)+y(33)+y(55));
residual(16) = lhs - rhs;
lhs = y(16);
rhs = (1-params(14))*(y(16)+y(74)-y(72));
residual(17) = lhs - rhs;
lhs = y(17);
rhs = (1-params(14))*(y(17)+y(75)-y(73));
residual(18) = lhs - rhs;
lhs = exp(y(45));
rhs = T(6)*T(12);
residual(19) = lhs - rhs;
lhs = exp(y(46));
rhs = params(37)*exp(y(76))+params(38)*exp(y(77));
residual(20) = lhs - rhs;
lhs = exp(y(76));
rhs = T(3)*T(13);
residual(21) = lhs - rhs;
lhs = exp(y(77));
rhs = T(3)*T(14);
residual(22) = lhs - rhs;
lhs = y(80);
rhs = params(1)*params(14)*(y(16)+y(80))+T(8)*(params(23)+y(34)-y(54));
residual(23) = lhs - rhs;
lhs = y(81);
rhs = params(1)*params(14)*(y(18)+y(81))+T(8)*(params(23)+y(54)*params(22)+y(34));
residual(24) = lhs - rhs;
lhs = y(16);
rhs = (1-params(14))*(y(16)+y(80)-y(78));
residual(25) = lhs - rhs;
lhs = y(18);
rhs = (1-params(14))*(y(18)+y(81)-y(79));
residual(26) = lhs - rhs;
lhs = exp(y(47));
rhs = params(39)*exp(y(82))+params(40)*exp(y(83))+params(41)*exp(y(84));
residual(27) = lhs - rhs;
lhs = exp(y(82));
rhs = T(10)*T(15);
residual(28) = lhs - rhs;
lhs = exp(y(83));
rhs = T(10)*T(16);
residual(29) = lhs - rhs;
lhs = exp(y(84));
rhs = T(10)*T(17);
residual(30) = lhs - rhs;
lhs = y(88);
rhs = params(1)*params(14)*(y(17)+y(88))+T(8)*(params(23)+y(55)+y(34)-y(54));
residual(31) = lhs - rhs;
lhs = y(89);
rhs = params(1)*params(14)*(y(18)+y(89))+T(8)*(params(23)+y(34)-params(22)*y(55));
residual(32) = lhs - rhs;
lhs = y(90);
rhs = params(1)*params(14)*(y(16)+y(90))+T(8)*(params(23)+y(34)-y(54)-params(22)*y(55));
residual(33) = lhs - rhs;
lhs = y(17);
rhs = (1-params(14))*(y(17)+y(88)-y(85));
residual(34) = lhs - rhs;
lhs = y(18);
rhs = (1-params(14))*(y(18)+y(89)-y(86));
residual(35) = lhs - rhs;
lhs = y(16);
rhs = (1-params(14))*(y(16)+y(90)-y(87));
residual(36) = lhs - rhs;
lhs = exp(y(48));
rhs = T(10)*T(18);
residual(37) = lhs - rhs;
lhs = exp(y(49));
rhs = params(42)*exp(y(91))+params(43)*exp(y(92));
residual(38) = lhs - rhs;
lhs = exp(y(91));
rhs = T(3)*T(19);
residual(39) = lhs - rhs;
lhs = exp(y(92));
rhs = T(3)*T(20);
residual(40) = lhs - rhs;
lhs = y(95);
rhs = params(1)*params(14)*(y(16)+y(95))+T(8)*(params(23)+y(35)-y(55));
residual(41) = lhs - rhs;
lhs = y(96);
rhs = params(1)*params(14)*(y(17)+y(96))+T(8)*(params(23)+params(22)*y(55)+y(35));
residual(42) = lhs - rhs;
lhs = y(16);
rhs = (1-params(14))*(y(16)+y(95)-y(93));
residual(43) = lhs - rhs;
lhs = y(17);
rhs = (1-params(14))*(y(17)+y(96)-y(94));
residual(44) = lhs - rhs;
lhs = exp(y(50));
rhs = params(44)*exp(y(97))+params(45)*exp(y(98))+params(46)*exp(y(99));
residual(45) = lhs - rhs;
lhs = exp(y(97));
rhs = T(6)*T(21);
residual(46) = lhs - rhs;
lhs = exp(y(98));
rhs = T(6)*T(22);
residual(47) = lhs - rhs;
lhs = exp(y(99));
rhs = T(6)*T(23);
residual(48) = lhs - rhs;
lhs = y(103);
rhs = params(1)*params(14)*(y(17)+y(103))+T(8)*(params(23)+y(35)-y(54)*params(22));
residual(49) = lhs - rhs;
lhs = y(104);
rhs = params(1)*params(14)*(y(18)+y(104))+T(8)*(params(23)+y(54)+y(35)-y(55));
residual(50) = lhs - rhs;
lhs = y(105);
rhs = params(1)*params(14)*(y(16)+y(105))+T(8)*(params(23)+y(35)-y(55)-y(54)*params(22));
residual(51) = lhs - rhs;
lhs = y(17);
rhs = (1-params(14))*(y(17)+y(103)-y(100));
residual(52) = lhs - rhs;
lhs = y(18);
rhs = (1-params(14))*(y(18)+y(104)-y(101));
residual(53) = lhs - rhs;
lhs = y(16);
rhs = (1-params(14))*(y(16)+y(105)-y(102));
residual(54) = lhs - rhs;
residual(55) = params(8)*y(1)+T(4)*(params(37)*y(78)+params(38)*(y(79)-y(54))+params(42)*y(93)+params(43)*(y(94)-y(55)));
residual(56) = params(8)*y(2)+T(4)*(params(33)*(y(66)+y(54))+params(34)*y(67)+params(46)*(y(54)+y(102))+params(44)*(y(54)+y(100)-y(55))+params(45)*y(101));
residual(57) = params(8)*y(3)+T(4)*(params(35)*(y(72)+y(55))+params(36)*y(73)+params(39)*y(85)+params(40)*(y(55)+y(86)-y(54))+params(41)*(y(55)+y(87)));
lhs = y(10);
rhs = params(1)*params(14)*(y(16)+y(10))+T(8)*(y(33)+params(23));
residual(58) = lhs - rhs;
lhs = y(11);
rhs = params(1)*params(14)*(y(18)+y(11))+T(8)*(params(23)+y(34));
residual(59) = lhs - rhs;
lhs = y(12);
rhs = params(1)*params(14)*(y(17)+y(12))+T(8)*(params(23)+y(35));
residual(60) = lhs - rhs;
lhs = y(16);
rhs = (1-params(14))*(y(16)+y(10)-y(1));
residual(61) = lhs - rhs;
lhs = y(18);
rhs = (1-params(14))*(y(18)+y(11)-y(2));
residual(62) = lhs - rhs;
lhs = y(17);
rhs = (1-params(14))*(y(17)+y(12)-y(3));
residual(63) = lhs - rhs;
lhs = y(5);
rhs = params(36)*y(73)+params(35)*y(72);
residual(64) = lhs - rhs;
lhs = y(9);
rhs = params(42)*y(93)+params(43)*y(94);
residual(65) = lhs - rhs;
lhs = y(4);
rhs = params(34)*y(67)+params(33)*y(66);
residual(66) = lhs - rhs;
lhs = y(6);
rhs = params(37)*y(78)+params(38)*y(79);
residual(67) = lhs - rhs;
lhs = y(7);
rhs = params(39)*y(85)+params(40)*y(86)+params(41)*y(87);
residual(68) = lhs - rhs;
lhs = y(8);
rhs = params(45)*y(101)+params(44)*y(100)+params(46)*y(102);
residual(69) = lhs - rhs;
lhs = y(51);
rhs = exp(y(42)+y(1))+params(35)*exp(y(70)+y(72))+params(36)*exp(y(71)+y(55)+y(73))+params(33)*exp(y(64)+y(66))+params(34)*exp(y(65)+y(67))-exp(y(39)+y(33));
residual(70) = lhs - rhs;
lhs = y(52);
rhs = exp(y(45)+y(2))+params(37)*exp(y(76)+y(54)+y(78))+params(38)*exp(y(77)+y(79))+params(41)*exp(y(84)+y(54)+y(87))+params(40)*exp(y(83)+y(86))+params(39)*exp(y(82)+y(54)+y(85)-y(55))-exp(y(41)+y(34));
residual(71) = lhs - rhs;
lhs = y(53);
rhs = exp(y(48)+y(3))+params(42)*exp(y(91)+y(55)+y(93))+params(43)*exp(y(92)+y(94))+params(44)*exp(y(97)+y(100))+params(46)*exp(y(99)+y(55)+y(102))+params(45)*exp(y(98)+y(55)+y(101)-y(54))-exp(y(40)+y(35));
residual(72) = lhs - rhs;
lhs = T(24);
rhs = T(24)*params(1)*(1+y(19))/(1+y(16));
residual(73) = lhs - rhs;
lhs = T(25);
rhs = T(25)*params(1)*(1+y(21))/(1+y(18));
residual(74) = lhs - rhs;
lhs = T(26);
rhs = T(26)*params(1)*(1+y(20))/(1+y(17));
residual(75) = lhs - rhs;
lhs = exp(y(15))+exp(y(54))*(1+y(22))*y(25)/(1+y(16));
rhs = y(52)+exp(y(54))*y(25)+exp(y(28)+y(31));
residual(76) = lhs - rhs;
lhs = exp(y(14))+exp(y(55))*(1+y(23))*y(26)/(1+y(16));
rhs = y(53)+exp(y(55))*y(26)+exp(y(29)+y(32));
residual(77) = lhs - rhs;
lhs = y(16);
rhs = params(1)*y(16)+T(27)*((y(13)-y(13)*params(32))*params(5)/(1-params(32))+params(6)*y(30)+log(params(7))-y(27)+params(13));
residual(78) = lhs - rhs;
lhs = y(18);
rhs = params(1)*y(18)+T(27)*(params(13)+log(params(7))+(y(15)-y(15)*params(32))*params(5)/(1-params(32))+y(31)*params(6)-y(28));
residual(79) = lhs - rhs;
lhs = y(17);
rhs = params(1)*y(17)+T(27)*(params(13)+log(params(7))+(y(14)-y(14)*params(32))*params(5)/(1-params(32))+y(32)*params(6)-y(29));
residual(80) = lhs - rhs;
lhs = y(33);
rhs = y(27)*(1-params(2))+log(1+y(19)*params(20))-y(59)-(params(2)*log(params(2))+(1-params(2))*log(1-params(2)));
residual(81) = lhs - rhs;
lhs = y(34);
rhs = y(28)*(1-params(2))+log(1+y(22)*params(20))-y(60)-(params(2)*log(params(2))+(1-params(2))*log(1-params(2)));
residual(82) = lhs - rhs;
lhs = y(35);
rhs = y(29)*(1-params(2))+log(1+y(23)*params(20))-y(61)-(params(2)*log(params(2))+(1-params(2))*log(1-params(2)));
residual(83) = lhs - rhs;
lhs = exp(y(39))*(1-params(2))/exp(y(30));
rhs = exp(y(27))/exp(y(33));
residual(84) = lhs - rhs;
lhs = exp(y(41))*(1-params(2))/exp(y(31));
rhs = exp(y(28))/exp(y(34));
residual(85) = lhs - rhs;
lhs = exp(y(40))*(1-params(2))/exp(y(32));
rhs = exp(y(29))/exp(y(35));
residual(86) = lhs - rhs;
lhs = exp(y(39))*params(2)/exp(y(36));
rhs = 1/exp(y(33));
residual(87) = lhs - rhs;
lhs = exp(y(41))*params(2)/exp(y(37));
rhs = 1/exp(y(34));
residual(88) = lhs - rhs;
lhs = exp(y(40))*params(2)/exp(y(38));
rhs = 1/exp(y(35));
residual(89) = lhs - rhs;
lhs = y(22);
rhs = y(19)+params(10)*(exp(y(25)-params(11))-1)+x(5);
residual(90) = lhs - rhs;
lhs = y(23);
rhs = y(19)+params(10)*(exp(y(26)-params(11))-1)+x(4);
residual(91) = lhs - rhs;
lhs = y(19)-params(9);
rhs = (y(19)-params(9))*params(53)+(1-params(53))*(y(16)*params(47)+params(50)*(y(39)-params(31)))+y(56);
residual(92) = lhs - rhs;
lhs = y(21)-params(9);
rhs = (y(21)-params(9))*params(54)+(1-params(54))*(y(18)*params(48)+params(51)*(y(41)-params(31)))+y(57);
residual(93) = lhs - rhs;
lhs = y(20)-params(9);
rhs = (y(20)-params(9))*params(17)+(1-params(17))*(y(17)*params(49)+params(52)*(y(40)-params(31)))+y(58);
residual(94) = lhs - rhs;
lhs = y(21)-y(22);
rhs = y(18)-y(16)+y(62);
residual(95) = lhs - rhs;
lhs = y(20)-y(23);
rhs = y(17)-y(16)+y(63);
residual(96) = lhs - rhs;
residual(97) = y(26)+y(25)+y(24);
lhs = y(56);
rhs = y(56)*params(18)+x(1);
residual(98) = lhs - rhs;
lhs = y(57);
rhs = y(57)*params(18)+x(2);
residual(99) = lhs - rhs;
lhs = y(58);
rhs = y(58)*params(18)+x(3);
residual(100) = lhs - rhs;
lhs = y(62);
rhs = y(62)*params(28)+x(9);
residual(101) = lhs - rhs;
lhs = y(63);
rhs = y(63)*params(28)+x(10);
residual(102) = lhs - rhs;
lhs = y(59)-params(25);
rhs = (y(59)-params(25))*params(26)+x(6);
residual(103) = lhs - rhs;
lhs = y(60)-params(25);
rhs = params(26)*(y(60)-params(25))+x(7);
residual(104) = lhs - rhs;
lhs = y(61)-params(25);
rhs = params(26)*(y(61)-params(25))+x(8);
residual(105) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
