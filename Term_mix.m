function t = Term_mix(a,ekh)
% This function calculates the summation terms as described in eq. 13 of
% main reference

% Load the function
load myFunc_mix.mat f
delta = exp(-ekh)/ekh;
w = ( 1 + ekh) * delta;
t = w * f(a,delta);
end