function [n,s] = TruncIndex(a,b,eh, kappa, nmax, s, eta)
% This function computes the initial guess

ekh = eh;
for n = 1 : nmax
    t = Term(a,b,ekh, eta);
    s = s + t;
    if abs(t) <= (kappa * abs(s))
        break;
    end
    ekh = ekh * eh;
end
end