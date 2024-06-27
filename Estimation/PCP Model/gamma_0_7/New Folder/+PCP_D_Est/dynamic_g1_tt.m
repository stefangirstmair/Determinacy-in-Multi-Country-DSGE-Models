function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 80);

T = PCP_D_Est.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(38) = (-params(21));
T(39) = params(8)*T(38)*getPowerDeriv(1-params(21)*y(39),T(1),1);
T(40) = params(8)*T(38)*getPowerDeriv(1-params(21)*y(40),T(1),1);
T(41) = params(8)*T(38)*getPowerDeriv(1-params(21)*y(41),T(1),1);
T(42) = (-(params(1)*params(14)));
T(43) = getPowerDeriv(exp(y(51)-params(32)*y(4)),(-params(5)),1);
T(44) = getPowerDeriv(exp(y(156)-y(51)*params(32)),(-params(5)),1);
T(45) = getPowerDeriv(exp(y(52)-params(32)*y(5)),(-params(5)),1);
T(46) = getPowerDeriv(exp(y(157)-y(52)*params(32)),(-params(5)),1);
T(47) = getPowerDeriv(exp(y(53)-params(32)*y(6)),(-params(5)),1);
T(48) = getPowerDeriv(exp(y(158)-y(53)*params(32)),(-params(5)),1);
T(49) = (-(exp(y(77))*(1-params(2))*exp(y(68))))/(exp(y(68))*exp(y(68)));
T(50) = (-(exp(y(79))*(1-params(2))*exp(y(69))))/(exp(y(69))*exp(y(69)));
T(51) = (-(exp(y(78))*(1-params(2))*exp(y(70))))/(exp(y(70))*exp(y(70)));
T(52) = (-((-(exp(y(65))*exp(y(71))))/(exp(y(71))*exp(y(71)))));
T(53) = (-((-(exp(y(66))*exp(y(72))))/(exp(y(72))*exp(y(72)))));
T(54) = (-((-(exp(y(67))*exp(y(73))))/(exp(y(73))*exp(y(73)))));
T(55) = (-(exp(y(74))*exp(y(77))*params(2)))/(exp(y(74))*exp(y(74)));
T(56) = (-(exp(y(75))*exp(y(79))*params(2)))/(exp(y(75))*exp(y(75)));
T(57) = (-(exp(y(76))*exp(y(78))*params(2)))/(exp(y(76))*exp(y(76)));
T(58) = T(4)*T(38)*getPowerDeriv(1-params(21)*(y(104)+y(92)),T(1),1);
T(59) = (-(T(6)*T(58)));
T(60) = getPowerDeriv(1-params(21)*(y(117)-y(92)),T(1),1);
T(61) = getPowerDeriv(1-params(21)*(y(93)+y(124)-y(92)),T(1),1);
T(62) = getPowerDeriv(1-params(21)*(y(92)+y(138)-y(93)),T(1),1);
T(63) = T(4)*T(38)*getPowerDeriv(1-params(21)*(y(92)+y(140)),T(1),1);
T(64) = (-(T(6)*T(63)));
T(65) = exp(y(92))*(1+y(10))*y(12)/(1+y(54))-exp(y(92))*y(63);
T(66) = T(4)*T(38)*getPowerDeriv(1-params(21)*(y(110)+y(93)),T(1),1);
T(67) = (-(T(10)*T(66)));
T(68) = T(4)*T(38)*getPowerDeriv(1-params(21)*(y(93)+y(125)),T(1),1);
T(69) = (-(T(10)*T(68)));
T(70) = getPowerDeriv(1-params(21)*(y(132)-y(93)),T(1),1);
T(71) = exp(y(93))*(1+y(11))*y(13)/(1+y(54))-exp(y(93))*y(64);
T(72) = (-(params(36)*exp(y(109)+y(93)+y(111))));
T(73) = T(4)*T(38)*getPowerDeriv(1-params(21)*y(105),T(1),1);
T(74) = T(4)*T(38)*getPowerDeriv(1-params(21)*y(111),T(1),1);
T(75) = T(4)*T(38)*getPowerDeriv(1-params(21)*y(116),T(1),1);
T(76) = (-(params(39)*exp(y(120)+y(92)+y(123)-y(93))));
T(77) = T(4)*T(38)*getPowerDeriv(1-params(21)*y(123),T(1),1);
T(78) = T(4)*T(38)*getPowerDeriv(1-params(21)*y(131),T(1),1);
T(79) = (-(params(45)*exp(y(136)+y(93)+y(139)-y(92))));
T(80) = T(4)*T(38)*getPowerDeriv(1-params(21)*y(139),T(1),1);

end
