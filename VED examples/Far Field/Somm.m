function y = Somm(c,d)
% This function generates the two-argument function.
% Specific to TM VED case only
%% Global Parameters
global i % index number of the distance array
global p % distance
global theta 
% Courtesy of Mazin M Mustafa

kp = c + d;
z_0 = 50; % height of the source
% Material Parameters
f = 10e6;
omega = 2*pi*f;
ep1 = 1;
% ep2 = 81 - 1i*7190.04;
ep2 = 4;
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);
eta1 = sqrt(mu0/ep0);

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

% TM X-tic impedance
Z1 = kz1/(omega*ep1);
Z2 = kz2/(omega*ep2);
Z3 = kz3/(omega*ep3);

% Stationary phase point
kp_s = k1 * sin(theta(i));
kz1_s = k1 * cos(theta(i));
k_s = sqrt(kp_s.^2 + kz1_s.^2);

pp = p * cos(theta(i));

gamma_1e = (kz2 / ep2 - kz1) / (kz2 / ep2 + kz1);

G_1 = (1 - gamma_1e) / (2i * kz1) * exp(-1i* kz1 * z_0); % Green's function for TM case

% Model f_theta from eq. 175 [1]
% for z directed source

f_theta =  - eta1^2 * mu0 /ep1 * sin(theta(i)) * sin(theta(i)) * G_1;

% the e-field from eq. 175
% considering a delta source along the z-axis

% E_theta = exp(1i * kz1 * z_0) * f_theta;

% TM case 
y = 4* pi* f_theta ;%* exp(1j* kz1_s * z_0) ;%* besselj(0, kp * p) * kp; % Sommerfeld Integrand for TM case

end
