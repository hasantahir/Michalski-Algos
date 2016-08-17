function [n,s] = TruncIndex_mix_PE(fun, a,eh, kappa, nmax, s)
% This function computes the initial guess
f = fun;
ekh = eh;
for n = 1 : nmax 
    t = Term_mix_PE(f, a,ekh);
    s = s + t;
    if abs(t) <= (kappa * abs(s))
        break;
    end
    ekh = ekh * eh;
end
end