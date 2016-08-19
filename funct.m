function y = funct(c,d)
% This function implements integrand of I_1(\rho)

% Courtesy of Mazin M Mustafa

x = c+d;

% When singularity lies ahead
if d >= 0
    y = ((1+x)^-0.5)*((1-x)^-0.5);
else
    y = ((1+x)^-0.5)*((-d)^-0.5);
end

% When singularity lies behind
%     if d <= 0
%         y = ((1+x)^-0.5)*((1-x)^-0.5);
%     end
%     if d > 0
%         y = (1/sqrt(-1))*((1+x)^-0.5)*((d)^-0.5);
%     end
% end
% p = pi/2;
% p = 53*pi/2;
% y = besselj(0,p*x)*x*y;
% Examples from
y = x/(1+x^2)*besselj(0,x);
y = x^2*besselj(0,x);
y = 1/2*log(x^2 + 1)*besselj(1,x);
y = (1-exp(-x))/(x*log(1 + sqrt(2)))*besselj(0,x);
end

% % Bessel Function reference
% @article{lucas1995evaluating,
%   title={Evaluating infinite integrals involving Bessel functions of arbitrary order},
%   author={Lucas, SK and Stone, HA},
%   journal={Journal of Computational and Applied Mathematics},
%   volume={64},
%   number={3},
%   pages={217--231},
%   year={1995},
%   publisher={Elsevier}
% }