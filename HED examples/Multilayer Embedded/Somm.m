function y = Somm(c,d)
% This function generates the two-argument function.
% Specific to HED case only
%% Global Parameters
global i % index number of the distance array
global p % distance
global nu % Switch for TE/TM case (alpha = 0 -> TE, else -> TM)

% Courtesy of Mazin M Mustafa

kp = c + d;

% Material Parameters
f = 1e12;
omega = 2*pi*f;
ep1 = 1; % Air
ep2 = 5.38; % GaN layer
ep3 = 11.7; % Silicon base
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);
k3 = omega*sqrt(mu0*ep0*ep3);

% Source location
zp = -50e-9 ;

% Observation location
z = -20e-9; % it has to be in the same layer as the source

% Layer Interface
z0 = -100e-9;
z1 = 0;
d = abs(z1 - z0);

kz1 = sqrt(k1 ^2 - kp ^2);

kz2 = sqrt(k2 ^2 - kp ^2);

kz3 = sqrt(k3 ^2 - kp ^2); 

% end



% Satisfying Radiation Condition
if imag(kz2) > 0
    kz2 = conj(kz2);
end
if imag(kz1) > 0
    kz1 = conj(kz1);
end
if imag(kz3) > 0
    kz3 = conj(kz3);
end

% Impedances
if nu == 0
    % TE case
    Z1 = omega/kz1;
    Z2 = omega/kz2;
    Z3 = omega/kz3;
else
    % TM case
    Z1 = kz1/(omega*ep1);
    Z2 = kz2/(omega*ep2);
    Z3 = kz3/(omega*ep3);
end


Gamma_left = (Z3 - Z2) / (Z3 + Z2);
Gamma_right = (Z1 - Z2) / (Z1 + Z2);

A = (Gamma_left * exp(1i*kz2*2*z0))/(1 - Gamma_left*Gamma_right*exp(-2i * kz2 * d)) ...
    *( exp(-1i * kz2 * zp) + Gamma_right * exp( -1i * kz2 * (2*z1 - zp)));

B = (Gamma_right * exp(-1i*kz2*2*z1))/(1 - Gamma_left*Gamma_right*exp(-2i * kz2 * d)) ...
    *( exp(+1i * kz2 * zp) + Gamma_left * exp( +1i * kz2 * (2*z0 - zp)));

% Define Voltage and Current based TL Greens Functions

V = Z2/2 * ( exp(-1i * kz2 * zp) + A*exp(-1i*kz2*z) + B*exp(+1i*kz2*z)); % TL GF 
if z < zp
    I = 1/2 * ( -exp(-1i * kz2 * zp) + A*exp(-1i*kz2*z) - B*exp(+1i*kz2*z)); % TL GF 
else
    I = 1/2 * ( exp(-1i * kz2 * zp) + A*exp(-1i*kz2*z) - B*exp(+1i*kz2*z)); % TL GF 
end


if nu == 0 % TE case
    y = V * besselj(nu, kp * p(i)) * kp; % Sommerfeld Integrand for TE case
else
    y = I * besselj(nu, kp * p(i)) * kp; % Green's function for TM case
    
end
% y = 1/(2i*kz1) * besselj(nu, kp * p(i)) * kp/(2*pi);

end
