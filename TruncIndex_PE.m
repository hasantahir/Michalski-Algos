function [n,s] = TruncIndex_PE(fun, a,b,eh, kappa, nmax, s, eta)
% This function computes the initial guess
f = fun;
ekh = eh;
for n = 1 : nmax 
    t = Term_PE(f, a,b,ekh, eta);
    s = s + t;
    if abs(t) <= (kappa * abs(s))
        break;
    end
    ekh = ekh * eh;
end
end