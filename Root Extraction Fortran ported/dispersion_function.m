function out = dispersion_function(gamma)
% ***********************************************************************
% 
%      Computes the dispersion function defined in the paper.
% 
% ***********************************************************************

% Data from datafile of test_case3 from the reference paper
lambda = 340; % in nanometers
k0  = 2*pi/lambda;

epsilon1 = 1.5;
epsilon2 = 1;
epsilonm = -0.291844 +1i*6.511223;

% Dimensions in nanometers
d = 50;
a1 = 50;
a2 = 50;

% Arguments of the upcoming trig functions
k0d  = k0*d;
k0a1 = k0*a1;
k0a2 = k0*a2;

gammasqd = gamma.^2;
% 
ztemp = gammasqd + epsilonm;

beta1 = sqrt(ztemp - epsilon1);
beta2 = sqrt(ztemp - epsilon2);

% 
% ---- Now compute the dispersion function for this gamma.
% 
%      Unnormalized form
% 
%      We start with the first two terms and their factor.
% 

factor12 = sinh(gamma*k0d)/gamma;

term1    = gammasqd .* cosh(k0a1 * beta1) .* cosh(k0a2 * beta2);
 
term2    = epsilonm*epsilonm.*beta1.*beta2/(epsilon1 * epsilon2);

term2    = term2 .* sinh(k0a1 * beta1) .* sinh(k0a2 * beta2);
% 
% .... Next, the third and fourth terms
% 
factor34 = cosh(k0d * gamma)*epsilonm;

term3    = beta2/epsilon2;
 
term3    = term3.*sinh(k0a2*beta2).*cosh(k0a1*beta1);

term4    = beta1/epsilon1;
 
term4    = term4.*sinh(k0a1*beta1).*cosh(k0a2*beta2);
% 
% .... Now put the whole lot together
% 
out = (factor12.*(term1 + term2)) + (factor34.*(term3 + term4));
end