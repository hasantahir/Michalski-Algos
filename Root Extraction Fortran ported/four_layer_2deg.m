%% solution for dispersion relation for 4 layer 2deg structure

clear all

f = .25e12; % Frequency

%% Constants
e0 = 8.85e-12; % Free-space permittivity
c = 3e8; % speed of light
e = 1.60218e-19; % electron charge
m_e = 9.109e-31; % Electron mass

% free space wavenumber
w = 2*pi*f;
k0 = 2*pi*f/c;

% layer material properties
e1 = 1; % Air
e2 = 11.7; % AlGaAs
e3 = 12.4; % GaAs

% Layer geometry
d1 = 20e-9;
d2 = 50e-9;

% Transverse propagation constants
kz1 = @(kx) sqrt(k0^2 - kx.^2);
kz2 = @(kx) sqrt(k0^2*e2 - kx.^2);
kz3 = @(kx) sqrt(k0^2*e3 - kx.^2);

% 
Z1 = @(kx) kx./(w*e1*e0);
Z2 = @(kx) kx./(w*e1*e0);
Z3 = @(kx) kx./(w*e1*e0);

Zup = @(kx) Z2(kx).*(Z1(kx) + 1i*Z2(kx).*tan(kx*d1))./ ...
    (Z2(kx) + 1i*Z1(kx).*tan(kx*d1));

Zdown = @(kx) 1i*Z2(kx).*tan(kx*d2);

f = @(kx) 1./Zup(kx) + 1./Zdown(kx);

r = newtzero(f,k0);
root = r
% neff = sqrt(ef - (r/k0).^2)

