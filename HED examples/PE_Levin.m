function val = PE_Levin(a, tol, q)

% q optimally should be pi/rho.
%% Initialize
global nu
kmax = 12;
X = zeros(1, kmax + 2);
A = zeros(1, kmax + 1);
B = zeros(1, kmax + 1);
% if nu == 0
%     X(1) = q(1);
% else
%     X(1) = q(2);
% end
X(1) = a;%q(1); % Corresponds to x(-1) in algo

q(1) = q(1) - a;
diff_q = diff(q);
% 
% % This line only for Sommerfeld Identity
% q = diff_q(1);
s = 0; % Initial sum to zero
% Initial guess
val = 1;


%% Main Alogirthm
% The partition-extrapolation algorithm PE calling Levin-Sidi and DE rule

for k = 2 : kmax + 2
    if nu == 0
        X(k) = X(k - 1) + diff_q(k+1);
    else
        X(k) = X(k - 1) + diff_q(k+1);
    end
    u = TanhSinhQuad(X(k - 1), X(k), tol);
    s = s + u;
    
    % Type of Levin Transformation
%     omega = u * (k - 2 + 1); % u- transformation
        omega = u;      % t- transformation
    %     if k > 1
    [val, A, B] = LevinSidi(k, s, omega, X, A, B);
    %     end
    if (k > 3 && abs(val - old) < tol * abs(val))
        break
    end
    old = val;
end

end