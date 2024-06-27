function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
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

assert(length(T) >= 38);

T = LCP_D.static_resid_tt(T, y, x, params);

T(28) = exp(y(13)-y(13)*params(32))*(1-params(32))*getPowerDeriv(exp(y(13)-y(13)*params(32)),(-params(5)),1);
T(29) = exp(y(14)-y(14)*params(32))*(1-params(32))*getPowerDeriv(exp(y(14)-y(14)*params(32)),(-params(5)),1);
T(30) = exp(y(15)-y(15)*params(32))*(1-params(32))*getPowerDeriv(exp(y(15)-y(15)*params(32)),(-params(5)),1);
T(31) = (-(T(6)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*(y(66)+y(54)),T(1),1)));
T(32) = getPowerDeriv(1-params(21)*(y(79)-y(54)),T(1),1);
T(33) = getPowerDeriv(1-params(21)*(y(55)+y(86)-y(54)),T(1),1);
T(34) = getPowerDeriv(1-params(21)*(y(54)+y(100)-y(55)),T(1),1);
T(35) = (-(T(6)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*(y(54)+y(102)),T(1),1)));
T(36) = (-(T(10)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*(y(72)+y(55)),T(1),1)));
T(37) = (-(T(10)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*(y(55)+y(87)),T(1),1)));
T(38) = getPowerDeriv(1-params(21)*(y(94)-y(55)),T(1),1);

end
