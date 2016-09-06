function y = Somm(c,d)
% This function generates the two-argument function.
% Specific to TM VED case only
%% Global Parameters
global i % index number of the distance array
global p % distance)

% Courtesy of Mazin M Mustafa

kp = c + d;
h = 5; % height of the source
% Material Parameters
f = 10e6;
omega = 2*pi*f;
ep1 = 1;
% ep2 = 4;
ep2 = 81 - 1i*7190.04;
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);


kz1 = sqrt(k1 ^2 - kp ^2);
kz2 = sqrt(k2 ^2 - kp ^2);

% end



% Satisfying Radiation Condition
if imag(kz2) > 0
    kz2 = conj(kz2);
end
if imag(kz1) > 0
    kz1 = conj(kz1);
end

gamma_1e = (kz2 / ep2 - kz1) / (kz2 / ep2 + kz1);

G_1 = (1 - gamma_1e) / (2i * kz1) * exp(-1i* kz1 * h); % Green's function for TM case

% TM case 
y = 4* pi* G_1 * besselj(0, kp * p(i)) * kp; % Sommerfeld Integrand for TM case

end
