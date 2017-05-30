function old = TanhSinhQuad(a, b, tol)
% This function is the DE rule designed to run on Horizontal Electric
% Dipole only
% For general purpose usage, go to the root folder
% Parameters
maxlev = 5;
eta = 1;
kappa = 1e-15;
nmax = 24;
global h

% Transformation to incorporate any limit other than 0-1
sigma = (b-a)/2;
gamma = (b+a)/2;


% This is actually hit-and-trial
eh = exp(h);



% First term g'(0)*f(gamma)
s = eta * Somm(gamma,0); % Call function to calculate first term
[n,s] = TruncIndex(a,b,eh, kappa, nmax, s, eta);
old = sigma * h * s;
for m  = 1 : maxlev
    e2h = eh;
    h = h/2;
    eh = exp(h);
    s =  PartSum(a,b,eh, e2h, n, eta);
    val = old/2 + sigma * h * s;
    if abs(val - old ) < (tol * abs(val))
%         disp('Converged');
        break;
    end
    old = val;
    n = 2 * n;
end
end
