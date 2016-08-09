function t = Term(a,b,ekh, eta)
% This function calculates the summation terms as described in eq. 13 of
% main reference

% Load the function
load myFunc.mat f
sigma = (b-a)/2;
q = exp(-eta * (ekh - 1/ekh));
delta = 2 * q / (1 + q);
w = eta *(ekh + 1/ekh) * delta/(1+q);   
t = w *( f(a, sigma*delta) + f(b, -sigma*delta) );
end