function t = Term(a,b,ekh, eta)
% This function generates the nodes and weights according to the tanh-sinh
% transform
% g(x) = sinh(x)
% t is the summation term 

% Load the function
sigma = (b-a)/2;
q = exp(-eta * (ekh - 1/ekh)); % q = exp( 2 *g(kh))
delta = 2 * q / (1 + q); % Expressing tanh
w = eta *(ekh + 1/ekh) * delta/(1+q);  % g'(t) *sech^2(g(t))
t = w *( funct(a, sigma*delta) + funct(b, -sigma*delta) ); %Term to summed
end