% This program implements Algorithm 7 from [1]
% The integrals should be of the type
% I_{\nv}(a, z, \rho) as defined in eq. (78)
clear all;close all
%% Global Parameters 
z = 1; % Second argument of the integral function
rho = 7.5; %Thirf argument of the integral function
q = pi/rho;
zeta = z;
nu = 0;
alpha = 1/2 - nu;
kmax = 10;
tol = 1e-15;

%% Load Function
load myFunc.mat f
%lower limit of the integral
a = 3.247; % Frist argument of the integral function
%% Initialize
X = zeros(1, kmax+2);
R = zeros(1, kmax+1);
mu = 2;
X(1) = a;
s = 0;

%% Main Alogirthm
for i = 1 : kmax + 1
    if i == 1
        X(i) = a + q;
        u = TanhSinhQuad(a, X(i), tol);
    else
        X(i) = X(i-1) + q;
        u = TanhSinhQuad(X(i-1), X(i), tol);
    end
    s = s + u;
    omega = Omega(i ,q, zeta , alpha, X);
    [val,R] = MosigMichalski(mu, i, s, omega, X, R);
    if (i > 1 && abs(val - old) < tol*abs(val))
        break
    end
    old = val;
end
val