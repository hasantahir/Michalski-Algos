function y = funct(c,d)
% This function implements integrand of I_1(\rho)

% Courtesy of Mazin M Mustafa

x = c+d;

    if d >= 0
        y = ((1+x)^-0.5)*((1-x)^-0.5);
    else
        y = ((1+x)^-0.5)*((-d)^-0.5);
    end
% else
%     if d <= 0
%         y = ((1+x)^-0.5)*((1-x)^-0.5);
%     end
%     if d > 0
%         y = (1/sqrt(-1))*((1+x)^-0.5)*((d)^-0.5);
%     end
% end
% p = pi/2;
p = 53*pi/2;
y = besselj(0,p*x)*x*y;
end