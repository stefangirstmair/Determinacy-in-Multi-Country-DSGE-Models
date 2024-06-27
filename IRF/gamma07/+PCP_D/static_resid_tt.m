function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 27);

T(1) = params(4)/params(21);
T(2) = params(8)*(1-params(21)*y(1))^T(1);
T(3) = exp(y(13))+exp(y(36));
T(4) = (1-params(8))/2;
T(5) = T(4)*(1-params(21)*(y(66)+y(54)))^T(1);
T(6) = exp(y(15))+exp(y(37));
T(7) = T(4)*(1-params(21)*y(67))^T(1);
T(8) = (1-params(1)*params(14))/(1+params(22));
T(9) = T(4)*(1-params(21)*(y(72)+y(55)))^T(1);
T(10) = exp(y(14))+exp(y(38));
T(11) = T(4)*(1-params(21)*y(73))^T(1);
T(12) = params(8)*(1-params(21)*y(2))^T(1);
T(13) = T(4)*(1-params(21)*y(78))^T(1);
T(14) = T(4)*(1-params(21)*(y(79)-y(54)))^T(1);
T(15) = T(4)*(1-params(21)*y(85))^T(1);
T(16) = T(4)*(1-params(21)*(y(55)+y(86)-y(54)))^T(1);
T(17) = T(4)*(1-params(21)*(y(55)+y(87)))^T(1);
T(18) = params(8)*(1-params(21)*y(3))^T(1);
T(19) = T(4)*(1-params(21)*y(93))^T(1);
T(20) = T(4)*(1-params(21)*(y(94)-y(55)))^T(1);
T(21) = T(4)*(1-params(21)*(y(54)+y(100)-y(55)))^T(1);
T(22) = T(4)*(1-params(21)*y(101))^T(1);
T(23) = T(4)*(1-params(21)*(y(54)+y(102)))^T(1);
T(24) = exp(y(13)-y(13)*params(32))^(-params(5));
T(25) = exp(y(15)-y(15)*params(32))^(-params(5));
T(26) = exp(y(14)-y(14)*params(32))^(-params(5));
T(27) = (1-params(1)*params(15))*(1-params(15))/(params(15)*(1+params(6)*params(12)));

end
