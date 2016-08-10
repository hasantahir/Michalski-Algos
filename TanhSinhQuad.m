function val = TanhSinhQuad(a, b, tol)
% This function 

% Parameters
eta = 1;
kappa = 1e-15;
nmax = 200;
maxlev = 5;

load myFunc.mat f

% 
sigma = (b-a)/2;
gamma = (b+a)/2;
% m = 0;
h = 1.5;
eh = exp(h);
s = eta * f(gamma,0);
[n,s] = TruncIndex(a,b,eh, kappa, nmax, s, eta);
old = sigma * h * s;
for m  = 1 : maxlev
    e2h = eh;
    h = h/2;
    eh = exp(h);
    s =  PartSum(a,b,eh, e2h, n, eta);
    val = old/2 + sigma * h * s;
    if abs(val - old ) < (tol * abs(val))
        break;
    end
    old = val;
    n = 2 * n;
end
end
