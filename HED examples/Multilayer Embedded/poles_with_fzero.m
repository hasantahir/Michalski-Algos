clear;close all;tic
% f = 1e12;
f = 1e12;
omega = 2*pi*f;
lambda = 3e8/f;
num = 1e2; %Size of the arrays

% Example Validations

% Material Properties
ep1 = 1; % Air
ep2 = 9 ; % GaN/AlGaN layers combined
ep3 = 11.9; % Silicon base

% EM constants
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

% Propagations Constants
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);
k3 = omega*sqrt(mu0*ep0*ep3);


% Middle Layer thickness
d = .5*lambda;

% Source Location
zp = -d/2;

% Layer Heights
z0 = -d;
z1 = 0;


% TE/TM switch
nu = 1;

kz1 = @(kp) sqrt(k1 ^2 - kp .^2);

kz2 = @(kp) sqrt(k2 ^2 - kp .^2);

kz3 = @(kp) sqrt(k3 ^2 - kp .^2);

% Define impedances
if nu == 0
    
    % TE Case
    Z1 = @(kp) omega./kz1(kp);
    Z2 = @(kp) omega./kz2(kp);
    Z3 = @(kp) omega./kz3(kp);
else
    % TM case
    Z1 = @(kp) kz1(kp)./(omega*ep1);
    Z2 = @(kp) kz2(kp)./(omega*ep2);
    Z3 = @(kp) kz3(kp)./(omega*ep3);
end

% Reflection Coefficients
Gamma_left = @(kp)(Z3(kp) - Z2(kp)) ./ (Z3(kp) + Z2(kp)); % Left-looking
Gamma_right = @(kp) (Z1(kp) - Z2(kp)) ./ (Z1(kp) + Z2(kp)); % Right-looking

% Unknown A
A = @(kp) (Gamma_left(kp) .* exp(1i*kz2(kp)*2*z0))./(1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * d)) ...
    .*( exp(-1i * kz2(kp) * zp) + Gamma_right(kp) .* exp( -1i * kz2(kp) * (2*z1 - zp)));

% Unknown B
B = @(kp) (Gamma_right(kp) .* exp(-1i*kz2(kp)*2*z1))./(1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * d)) ...
    .*( exp(+1i * kz2(kp) * zp) + Gamma_left(kp) .* exp( +1i * kz2(kp) * (2*z0 - zp)));

% Denominator
D = @(kp) 1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * d);
D1 = @(kp) D(kp)/k1;
x0 = [k1/k1 k3/k1];
z = fzero(Gamma_left,k2)