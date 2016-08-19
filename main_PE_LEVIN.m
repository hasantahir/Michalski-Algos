% This program implements Algorithm 7 from [1]
% The integrals should be of the type
% I_{\nv}(a, z, \rho) as defined in eq. (78)
clear all;close all
%% Global Parameters
z = 1; % Second argument of the integral function
rho = 1; % 1 seems to be the optimal value
q = pi/rho; % Discrete increment
nu = 0;
tol = 1e-15;

% lower limit of the integral
a = 0;

% Call PE routine
val = PE_Levin(a, tol, q)

