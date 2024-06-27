function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 31);

T(1) = params(4)/params(21);
T(2) = params(8)*(1-params(21)*y(39))^T(1);
T(3) = exp(y(51))+exp(y(74));
T(4) = (1-params(8))/2;
T(5) = T(4)*(1-params(21)*(y(104)+y(92)))^T(1);
T(6) = exp(y(53))+exp(y(75));
T(7) = T(4)*(1-params(21)*y(105))^T(1);
T(8) = (1-params(1)*params(14))/(1+params(22));
T(9) = T(4)*(1-params(21)*(y(110)+y(93)))^T(1);
T(10) = exp(y(52))+exp(y(76));
T(11) = T(4)*(1-params(21)*y(111))^T(1);
T(12) = params(8)*(1-params(21)*y(40))^T(1);
T(13) = T(4)*(1-params(21)*y(116))^T(1);
T(14) = T(4)*(1-params(21)*(y(117)-y(92)))^T(1);
T(15) = T(4)*(1-params(21)*y(123))^T(1);
T(16) = T(4)*(1-params(21)*(y(93)+y(124)-y(92)))^T(1);
T(17) = T(4)*(1-params(21)*(y(93)+y(125)))^T(1);
T(18) = params(8)*(1-params(21)*y(41))^T(1);
T(19) = T(4)*(1-params(21)*y(131))^T(1);
T(20) = T(4)*(1-params(21)*(y(132)-y(93)))^T(1);
T(21) = T(4)*(1-params(21)*(y(92)+y(138)-y(93)))^T(1);
T(22) = T(4)*(1-params(21)*y(139))^T(1);
T(23) = T(4)*(1-params(21)*(y(92)+y(140)))^T(1);
T(24) = exp(y(147)-y(51)*params(32))^(-params(5));
T(25) = params(1)*(1+y(57))*T(24);
T(26) = exp(y(149)-y(53)*params(32))^(-params(5));
T(27) = params(1)*(1+y(59))*T(26);
T(28) = exp(y(148)-y(52)*params(32))^(-params(5));
T(29) = params(1)*(1+y(58))*T(28);
T(30) = (1-params(1)*params(15))*(1-params(15))/(params(15)*(1+params(6)*params(12)));
T(31) = params(5)/(1-params(32));

end
