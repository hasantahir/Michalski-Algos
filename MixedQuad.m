function val = MixedQuad(a, tol)
% This function 

% Parameters
kappa = 1e-15;
nmax = 24;
maxlev = 5;

load myFunc_mix.mat f

% 
m = 0;
h = 1;
eh = exp(h);
delta = exp(-1);
w = 2 * delta;
s = w * f(a,delta);
[n1,s] = TruncIndex_mix(a,eh, kappa, nmax, s);
[n2,s] = TruncIndex_mix(a,1/eh, kappa, nmax, s);
old = h * s;
for m  = 1 : maxlev
    e2h = eh;
    h = h/2;
    eh = exp(h);
    s1 =  PartSum_mix(a,eh, e2h, n1);
    s2 =  PartSum_mix(a,1/eh, 1/e2h, n2);
    val = old/2 + h * (s1 + s2);
    if abs(val - old ) < (tol * abs(val))
        break;
    end
    old = val;
    n1 = 2 * n1;
    n2 = 2 * n2;
end
end
