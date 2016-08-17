function t = Term_mix(a,ekh)
% This function calculates the summation terms as described in eq. 13 of
% main reference

% Load the function
delta = exp(-ekh)/ekh;
w = ( 1 + ekh) * delta;
t = w * funct_mix(a,delta);

end