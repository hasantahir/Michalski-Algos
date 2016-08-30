function a = FindFirstZero()

% This function returns the first zero of the integrand function. This
% returned value will be set as the breakpoint to shift from one algo to
% another.

% The function defined in funct.m should be the same as here. However,
% we'll regenerate it as a function handle since the routine AllZeros takes
% only function handles.

% The first zero is chosen.
global p i z
% For figure 7b
t = 1;


% Lucas Decomposition
if p(i) == 0
    a = 0;
else
    Y_p = @(x) bessely(0,p(i)*x);
    Y_v = @(x) bessely(3/2,x*t);
    zero_yp = AllZeros(Y_p, 0, 20); % Find zeros of Y_
    zero_yv = AllZeros(Y_v, 0, 20);
    a = min(zero_yp(1)/p(i), zero_yv(1)/t);
%     a = zros(1);
end

end