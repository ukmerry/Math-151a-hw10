% qSplineEval:
% 
% A program for evaluating a quadradic spline that 
% has been constructed using n+1 equispace knots 
% over an interval [xMin,xMax]
%
% Inputs:
%
%       x  : vector of values at which the spline
%            is to be evaluated
%
% a,b,c    : vectors of spline coefficients assuming a spline
%            representation over the ith panel as
%
%           S_i(x) = a_i + b_i*(x-x_(i-1)) + c_i*(x-x_(i-1))^2
%          
% xMin,xMax : interval over which the spline was constructed 
%
% Output:
%
%  y        : vector of values of the spline evaluated at locations
%             specified by x  
% 
% Math 151A, Winter 2018 (03/08/2018)
% 
function [y] = qSplineEval(x,a,b,c,xMin,xMax)
  n = length(a);
  h = (xMax-xMin)/n;
  
  for i = 1:length(x)
  panelIndex = floor((x(i) - xMin)/h) + 1;
  if(panelIndex < 1) panelIndex = 1; end
  if(panelIndex > n) panelIndex = n; end
  
  xBase   = xMin + (panelIndex-1)*h;
  y(i)    = a(panelIndex) + b(panelIndex)*(x(i)- xBase) + c(panelIndex)*(x(i)-xBase)^2;
  end