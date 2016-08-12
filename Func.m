function Func()
% The integrand is defined as a function handle

% To avoid end-point singularity, the function should have two arguments

% Example 1 - taken from 
% http://www.codeproject.com/Articles/31550/Fast-Numerical-Integration
% f = @(c,d) exp(-(c+d)/5.0).*(2.0 + sin(2.0*(c+d)));

% Example 2 - taken from reference paper
% p = 53*pi/2;
% f = @(c,d) besselj(0,p*(c+d)).*(c+d)./sqrt(1 - (c+d).^2);

% % Example 3 - taken from the paper
% p = .011;
% f = @(c,d) exp(-sqrt((c+d).^2-1)*p)./sqrt((c+d).^2 -1).*(c+d);


% Example 4 - taken from http://keisan.casio.com/exec/system/1285343955
% f =@(c,d) sqrt((c+d).^2 - 1);

% Example 5 - Numerical Recipes pg 217
p = 1;
f = @(c,d) besselj(0,p*(c+d)).*(c+d)./(1 + (c+d).^2);
% f = @(x) besselj(0,(x)).*(x)./(1 + (x).^2);

save myFunc.mat f
save myFunc_mix.mat f
% To reuse it in other functions, load the mat file there
end