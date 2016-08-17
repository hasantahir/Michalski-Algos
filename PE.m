function val = PE(fun, a, tol, q, z, alpha)

%% Initialize
kmax = 10;
X = zeros(1, kmax + 2);
R = zeros(1, kmax + 1);
mu = 2;
X(1) = a;
s = 0;

%% Load Function
f = fun;
%% Main Alogirthm
for k = 2 : kmax + 2
    
    X(k) = X(k-1) + q;
    u = TanhSinhQuad_PE(f, X(k-1), X(k), tol);
    s = s + u;
    omega = Omega(k, q, z, alpha, X);
    [val, R] = MosigMichalski(mu, k, s, omega, X, R);
    if (k > 2 && abs(val - old) < tol*abs(val))
        break
    end
    old = val;
end
end