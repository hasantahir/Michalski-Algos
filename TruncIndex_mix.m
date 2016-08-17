function [n,s] = TruncIndex_mix(a,eh, kappa, nmax, s)
% This function computes the initial guess

ekh = eh;
for n = 1 : nmax 
    t = Term_mix(a,ekh);
    s = s + t;
    if abs(t) <= (kappa * abs(s))
        disp('Converged');
        break;
    end
    ekh = ekh * eh;
end
end