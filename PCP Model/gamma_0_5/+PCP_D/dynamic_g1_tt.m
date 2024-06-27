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

assert(length(T) >= 49);

T = PCP_D.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(36) = getPowerDeriv(exp(y(51)-params(32)*y(4)),(-params(5)),1);
T(37) = getPowerDeriv(exp(y(147)-y(51)*params(32)),(-params(5)),1);
T(38) = getPowerDeriv(exp(y(52)-params(32)*y(5)),(-params(5)),1);
T(39) = getPowerDeriv(exp(y(148)-y(52)*params(32)),(-params(5)),1);
T(40) = getPowerDeriv(exp(y(53)-params(32)*y(6)),(-params(5)),1);
T(41) = getPowerDeriv(exp(y(149)-y(53)*params(32)),(-params(5)),1);
T(42) = (-(T(6)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*(y(104)+y(92)),T(1),1)));
T(43) = getPowerDeriv(1-params(21)*(y(117)-y(92)),T(1),1);
T(44) = getPowerDeriv(1-params(21)*(y(93)+y(124)-y(92)),T(1),1);
T(45) = getPowerDeriv(1-params(21)*(y(92)+y(138)-y(93)),T(1),1);
T(46) = (-(T(6)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*(y(92)+y(140)),T(1),1)));
T(47) = (-(T(10)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*(y(110)+y(93)),T(1),1)));
T(48) = (-(T(10)*T(4)*(-params(21))*getPowerDeriv(1-params(21)*(y(93)+y(125)),T(1),1)));
T(49) = getPowerDeriv(1-params(21)*(y(132)-y(93)),T(1),1);

end
