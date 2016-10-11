% This program computes the Sommerfeld integral for a Vertical Electric
% Dipole
clear all; close all
tic
tol = 1e-6; % tolerance of the routine
num = 60; %Size of the arrays
%% Global Parameters
global i % index number of the distance array
global p % distance
global a % Breakpoint location
global maxlev %Settings for DE integration routine
global theta % Theta scan angle

f = 10e6;
c = 3e8;
lambda = c/f;
omega = 2*pi*f;

ep1 = 1;
% ep2 = 81 - 1i*7190.04;
ep2 = 4;
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);

a = 2*k1; % Set breakpoint
p = lambda * 1e2; % Define distance array
q = pi/p;

% p = 1e3/lambda;
val_1 = zeros(size(p));
val_2 = zeros(size(p));
val_3 = zeros(size(p));
val = zeros(size(p));

% This program computes the far-field of a vertical electric dipole
phi = linspace(0,2*pi,360);
theta = linspace( 0 + pi/70 , pi - pi/70 , num);

% for theta sweep
kp_s = k1 *sin(theta); % Stationary phase point
kzn_s = k1 * cos(theta); % Stationary phase point

% Integration routine
for i = 1 : length(theta)
        maxlev = 25;
        val_1(i) = TanhSinhQuad(0, k1 + .1i/p, tol); % Integrate upto k through DE
        val_2(i) = TanhSinhQuad(k1 + .1i/p, a, tol); % Integrate k upto a through DE
        maxlev = 15;
        val_3(i) = PE_Levin(a, tol, q); % Tail through PE Levin with Lucas
    
    val(i) = val_1(i) + val_2(i) + val_3(i);
end

% Radiation pattern
polar(theta, abs(val)./max(max(abs(val))))