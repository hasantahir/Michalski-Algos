function y = funct(c,d)
% This function implements integrand of I_1(\rho)

%% For sweep, uncomment these global variables
global z i p l

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
% % Examples from
% y = x/(1+x^2)*besselj(0,x);
% y = x^2*besselj(0,x);
% y = 1/2*log(x^2 + 1)*besselj(1,x);
% y = (1-exp(-x))/(x*log(1 + sqrt(2)))*besselj(0,x);

% from reference paper [1] eq. 80

% For figure 7a
% p = 0;
% t = 1;
% y = exp(-x*z(i))*besselj(0,p*x)*besselj(3/2,x*t);

% %% Reuse this part in FindFirstZeros.m
% % For figure 7b
% t = 1;
% % Lucas Decomposition
% J_plus = 1/2*(besselj(0,p(i)*x)*besselj(3/2,x*t) - ...
%     bessely(0,p(i)*x)*bessely(3/2,x*t));
% J_minus = 1/2*(besselj(0,p(i)*x)*besselj(3/2,x*t) + ...
%     bessely(0,p(i)*x)*bessely(3/2,x*t));
% J = besselj(0,p(i)*x)*besselj(3/2,x*t);
% if l == 1
%     y = exp(-x*z)*J_plus*x^.5;
% elseif l == 2
%     y = exp(-x*z)*J_minus*x^.5;
% else
%     y = exp(-x*z)*J*x^.5;
%     
% end


%% From reference [1], eq 78
% I_2(5.13562, 0, 1)
% z = 0;
% p = 1;
% a = 5.13562;
% y = besselj(2,x*p)*x^2;
% % y = @(x) besselj(2,p*x).*x.^2;


%% From reference [1], eq 78
% % I_0(3.247, 0, 7.5)
% z = 0;
% p = 7.5;
% % a = 3.247;
% a = .3206;
% y = besselj(0,x*p)*x^0;
% % % y = @(x) besselj(0,p*x).*x.^0;

%% From reference [1], eq 78
% I_0(3.247, 0, 7.5)
z = 1;
p = 7.5;
a = 3.247;
y = exp(-x*z)*besselj(0,x*p)*x^0;
% % y = @(x) exp(-x*z).*besselj(0,p*x).*x.^0;

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