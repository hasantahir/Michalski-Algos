function old = TanhSinhQuad(a, b, tol)
% This function 

% Parameters
eta = 1;
kappa = 1e-15;
nmax = 193;
maxlev = 5;


% Transformation to incorporate any limit other than 0-1
sigma = (b-a)/2;
gamma = (b+a)/2;


h = 1.5; % This is actually hit-and-trial
eh = exp(h);

% First term g'(0)*f(gamma)
s = eta * funct(gamma,0); % Call function to calculate first term
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
