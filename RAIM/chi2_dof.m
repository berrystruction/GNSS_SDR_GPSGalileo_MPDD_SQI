%                             chi2_dof.m
%  Scope:   This MATLAB macro computes the probability density function of the 
%           chi-square distribution with dof degrees of freedom.
%  Usage:   xoutput = chi2_dof(x,y)
%  Description of parameters:
%           x       - input, scalar independent variable
%           dof     - input, scalar, degree of freedom, global variable
%           xoutput - output, computed probability density function 
%  Remark:
%           The computed function is equivalent to the state derivative of the 
%           chi-square distribution with dof degrees of freedom.
%  Last update:  04/07/00
%  Copyright (C) 1997-00 by LL Consulting. All Rights Reserved.

function  xoutput = chi2_dof(x,dof)

k2 = 0.5*dof;
gfnc = gamma(k2);
power2k = 2^k2;

xoutput = x^(k2 - 1) * exp(-0.5*x) / power2k / gfnc;

end

