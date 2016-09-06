function val = PE_Levin(a, tol, q)

% q optimally should be pi/rho.
%% Initialize
kmax = 51;
X = zeros(1, kmax + 2);
A = zeros(1, kmax + 1);
B = zeros(1, kmax + 1);
X(1) = a; % Corresponds to x(-1) in algo
s = 0; % Initial sum to zero
% Initial guess
val = 1;

%% Main Alogirthm
% The partition-extrapolation algorithm PE calling Levin-Sidi and DE rule

for k = 2 : kmax + 2
    
    X(k) = X(k - 1) + q;
    u = TanhSinhQuad(X(k - 1), X(k), tol);
    s = s + u;
    
    % Type of Levin Transformation
<<<<<<< HEAD
%     omega = u * (k - 2 + 1); % u- transformation
        omega = u;      % t- transformation
=======
    omega = u * (k - 2 + 1); % u- transformation
    %     omega = u;      % t- transformation
>>>>>>> master
    %     if k > 1
    [val, A, B] = LevinSidi(k, s, omega, X, A, B);
    %     end
    if (k > 3 && abs(val - old) < tol * abs(val))
        break
    end
    old = val;
end

end