function val = PE_Levin(fun, a, tol, q)

%% Initialize
kmax = 10;
X = zeros(1, kmax + 2);
A = zeros(1, kmax + 1);
B = zeros(1, kmax + 1);
X(1) = a;
s = 0;

%% Load Function
f = fun;
%% Main Alogirthm
for k = 2 : kmax + 2
    X(k) = X(k-1) + q;
    u = TanhSinhQuad_PE(f, X(k-1), X(k), tol);
    s = s + u;
    omega = u;
    [val, A, B] = LevinSidi(k, s, omega, X, A, B);
    if (k > 2 && abs(val - old) < tol*abs(val))
        break
    end
    old = val;
end

end