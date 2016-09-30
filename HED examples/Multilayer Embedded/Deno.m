function D = Deno(kp)

% F(1) = exp(-exp(-(x(1)+x(2)))) - x(2)*(1+x(1)^2);
% F(2) = x(1)*cos(x(2)) + x(2)*sin(x(1)) - 0.5;
global nu

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

kz1 = sqrt(k1 ^2 - kp ^2);

kz2 = sqrt(k2 ^2 - kp ^2);

kz3 = sqrt(k3 ^2 - kp ^2);

% Define impedances
if nu == 0
    
    % TE Case
    Z1 = omega/kz1;
    Z2 = omega/kz2;
    Z3 = omega/kz3;
else
    % TM case
    Z1 = kz1/(omega*ep1);
    Z2 = kz2/(omega*ep2);
    Z3 = kz3/(omega*ep3);
end

% Reflection Coefficients
Gamma_left = (Z3 - Z2) / (Z3 + Z2); % Left-looking
Gamma_right = (Z1 - Z2) / (Z1 + Z2); % Right-looking

% Unknown A
A = (Gamma_left * exp(1i*kz2*2*z0))/(1 - Gamma_left*Gamma_right*exp(-2i * kz2 * d)) ...
    .*( exp(-1i * kz2 * zp) + Gamma_right .* exp( -1i * kz2 * (2*z1 - zp)));

% Unknown B
B =  (Gamma_right * exp(-1i*kz2*2*z1))/(1 - Gamma_left*Gamma_right*exp(-2i * kz2 * d)) ...
    *( exp(+1i * kz2 * zp) + Gamma_left * exp( +1i * kz2 * (2*z0 - zp)));

% Denominator
D =  1 - Gamma_left*Gamma_right*exp(-2i * kz2 * d);



end