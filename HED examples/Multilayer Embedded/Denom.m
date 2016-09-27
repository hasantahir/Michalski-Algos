function D = Denom(p)
% This function generates the denominator to find poles of the HED embedded
% in multilayer envorinoment
%% Global Parameters

global nu % Switch for TE/TM case (alpha = 0 -> TE, else -> TM)
global f
global ep1 ep2 ep3 
global k1 k2 k3
global mu0 ep0 
% Courtesy of Mazin M Mustafa

kp = p;

% Material Parameters

omega = 2*pi*f;
lambda = 3e8/f;


% Middle Layer thickness
% d = 1.02*lambda;
d = .8e-6;
% d = 0;

% Source Location
zp = -d/2;
% zp = 0;

% Layer Heights
z0 = -d;
z1 = 0;


kz1 = sqrt(k1 ^2 - kp ^2);

kz2 = sqrt(k2 ^2 - kp ^2);

kz3 = sqrt(k3 ^2 - kp ^2);

% Satisfying Radiation Condition
if imag(kz1) > 0
    kz1 = conj(kz1);
end

if imag(kz2) > 0
    kz2 = conj(kz2);
end

if imag(kz3) > 0
    kz3 = conj(kz3);
end

% Define impedances
if nu == 0
    
    % TE Case
    Z1 = omega*mu0/kz1;
    Z2 = omega*mu0/kz2;
    Z3 = omega*mu0/kz3;
else
    % TM case
    Z1 = kz1/(omega*ep1*ep0);
    Z2 = kz2/(omega*ep2*ep0);
    Z3 = kz3/(omega*ep3*ep0);
end

% Reflection Coefficients
Gamma_left = (Z3 - Z2) / (Z3 + Z2); % Left-looking
% For a grounded case
% Gamma_12 = 
Gamma_right =  (Z1 - Z2) / (Z1 + Z2); % Right-looking

% Unknown A
A =  (Gamma_left * exp(1i * kz2 * 2 * z0)) / (1 - Gamma_left * Gamma_right * exp(-2i * kz2 * d)) ...
    * ( exp(-1i * kz2 * zp) + Gamma_right * exp( -1i * kz2 * (2*z1 - zp)));

% Unknown B 
B = (Gamma_right * exp(-1i * kz2 * 2 * z1) ) / ( 1 - Gamma_left * Gamma_right * exp(-2i * kz2 * d)) ...
     * ( exp(+1i * kz2 * zp) + Gamma_left * exp( +1i * kz2 * (2 * z0 - zp)));

% Denominator
D = 1 - Gamma_left * Gamma_right * exp(-2i * kz2 * d);

end