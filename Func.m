function Func()
% The integrand is defined as a function handle

% To avoid end-point singularity, the function should have two arguments

f = @(c,d) exp(-(c+d)/5.0).*(2.0 + sin(2.0*(c+d)));

save myFunc.mat f
% To reuse it in other functions, load the mat file there
end