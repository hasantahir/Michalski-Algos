function y = funct_mix(c,d)
% This function implements integrand of I_2(\rho)

% Courtesy of Mazin M Mustafa

x = c + d;


% y = ((1+x)^-0.5)*((d)^-0.5);
% 
% z = .011;
% y = exp(-sqrt(x^2-1)*z)*x*y;
% From Online reference on bessel functions
y = x^2*besselj(0,x);
% y = x/(1+x^2)*besselj(100,x);

end