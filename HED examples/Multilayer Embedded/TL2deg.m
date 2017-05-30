function D = TL2deg(kp)
% ***********************************************************************
%
%      Computes the dispersion relation of a backgated HEMT with 
%      2DEG modeled as a charged sheet
%      
% ***********************************************************************


%% Calls polee_search_MIM
% Data from datafile of test_case3 from the reference paper
f = 25e12;
omega = 2*pi*f;
lambda = 3e8/f;
num = 80; %Size of the arrays

% Example Validations

% Material Properties
% ep1 = 1; % Air
% ep2 = 9.7 ; % GaN/AlGaN layers combined
% ep3 = 1; % Silicon base
% % backed by a PEC layer

d = 20e-9;
h = 50e-9;
% Chinese Homotopy Method
ep1 = 1; % Air
ep2 = 9.6 ; % AlGaAs layers combined
ep3 = 10; % GaAS
% EM constants
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

% Propagations Constants
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);
k3 = omega*sqrt(mu0*ep0*ep3);

%% Define 2DEG
q = 1.6012e-19; % Free elecrton charge
m = 9.109e-31; % Free electron mass

N = 2e11*1e4; % Electon Density % Sample D
mu = 15.6*1e6*1e-4; % Mobility % Sample D 
ms = .063*m; % Effective mass
tau = mu*ms/q; % Average Scattering time
% Conductivity
cond = N*q^2*tau/ms./(1 - 1i*omega*tau);
% 

% TE/TM switch

kz1 = (1*sqrtbr(k1 ^2 - kp .^2, pi/2));

kz2 = (+sqrt(k2 ^2 - kp .^2));

kz3 = sqrt(k3 ^2 - kp .^2);

   
    Z1 =  kz1./(omega*ep1);
    Z2 =  kz2./(omega*ep2);
    Z3 =  kz3./(omega*ep3);


% Normalize Impedances to free-space
% Z0 = sqrt(mu0/ep0);
% Z1 =  Z1/Z0;
% Z2 =  Z2/Z0;
% Z3 =  Z3/Z0;

%% Reflection Coefficients
Z_left = Z2 .* (Z1 + 1i*Z2.*tan(kz2 *d))./(Z2 + 1i*Z1.*tan(kz2 *d));
Z_right =  Z3 * 1i.*tan(kz3*h); % Right-looking

% Denominator, Dispersion relation
D =  1./Z_left + 1./Z_right - cond;

end