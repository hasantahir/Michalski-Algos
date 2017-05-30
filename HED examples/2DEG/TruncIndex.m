function [n,s] = TruncIndex(a,b,eh, kappa, nmax, s, eta)
% Actual answer of summation plus the first term is stored here

ekh = eh; % first term of summation, k = 1
for n = 1 : nmax 
    t = Term(a,b,ekh, eta); % Summation
    s = s + t; % Add first term and the summation
    if abs(t) <= (kappa * abs(s))
%         disp('Converged');
        break;
    end
    ekh = ekh * eh; %k*h = h*(k+1) 
end 
end