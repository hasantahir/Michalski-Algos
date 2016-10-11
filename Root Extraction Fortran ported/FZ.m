function D = FZ(kp)
  % !***********************************************************************
  % !
  % !     Computes the function "f" at the point "z".
  % !
  % !***********************************************************************
  % !
  % !---- Tell the compiler to use the module which computes the dispersion
  % !     function.
  % !
  % use dispersion_function_module
  % !
  % !---- All variables in this unit have explicit type
  % !

%     D = 54.0 + kp.*(44.0 + kp.*(20.0 - kp.*(3.0-kp)));
%     D = kp;
%     D = log - kp + 1.195281 + 0.536353*1i;
%   D = kp.^(-0.5*kp) - kp + 2.195093245 -2.750498700*1i;
% D = kp.^50 + kp.^12 - 5*sin(20*kp).*cos(12*kp) - 1;
% D = kp.^2 - 1;

% f = 1e12;
% omega = 2*pi*f;
% lambda = 3e8/f;
% 
% % Example Validations
% 
% % Material Properties
% ep1 = 1; % Air
% ep2 = 7.5 ; % GaN/AlGaN layers combined
% ep3 = 11; % Silicon base
% 
% % EM constants
% mu0 = 4*pi*1e-7;
% ep0 = 8.854e-12;
% 
% % Propagations Constants
% k1 = omega*sqrt(mu0*ep0*ep1);
% k2 = omega*sqrt(mu0*ep0*ep2);
% k3 = omega*sqrt(mu0*ep0*ep3);
% 
% 
% % Middle Layer thickness
% d = 1.5*lambda;
% 
% % Source Location
% zp = -d/2;
% 
% % Layer Heights
% z0 = -d;
% z1 = 0;
% 
% 
% % TE/TM switch
% nu = 1;
% 
% kz1 = sqrt(k1 ^2 - kp .^2);
% 
% kz2 = sqrt(k2 ^2 - kp .^2);
% 
% kz3 = sqrt(k3 ^2 - kp .^2);
% 
% % Define impedances
% if nu == 0
%     
%     % TE Case
%     Z1 = omega./kz1;
%     Z2 = omega./kz2;
%     Z3 = omega./kz3;
% else
%     % TM case
%     Z1 =  kz1./(omega*ep1);
%     Z2 =  kz2./(omega*ep2);
%     Z3 = kz3./(omega*ep3);
% end
% 
% Left Looking Input Impedance
% Z_in_left = Z2 * (Z3 + Z2*1i*tan(kz2*d/2))/(Z2 + Z3*1i*tan(kz2*d/2));
% 
% % Right Looking Input Impedance
% Z_in_right = Z2 * (Z1 + Z2*1i*tan(kz2*d/2))/(Z2 + Z1*1i*tan(kz2*d/2));

% % Reflection Coefficients
% Gamma_left = (Z_in_left - Z2) ./ (Z_in_left + Z2); % Left-looking
% Gamma_right =  (Z_in_right - Z2) ./ (Z_in_right + Z2); % Right-looking
% 
% % Unknown A
% A =  (Gamma_left .* exp(1i*kz2*2*z0))./(1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * d)) ...
%     .*( exp(-1i * kz2 * zp) + Gamma_right .* exp( -1i * kz2 * (2*z1 - zp)));
% 
% % Unknown B
% B = (Gamma_right .* exp(-1i*kz2*2*z1))./(1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * d)) ...
%     .*( exp(+1i * kz2 * zp) + Gamma_left .* exp( +1i * kz2 * (2*z0 - zp)));
% 
% % Denominator
% D =  1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * d);

%% Example from Chinese paper
f = 5e9;
omega = 2*pi*f;
lambda = 3e8/f;
num = 5e2; %Size of the arrays

% Example Validations

% Material Properties
ep1 = 12; 
ep2 = 4 ; 
ep3 = 2.1;
ep4 = 1;


% EM constants
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

% Propagations Constants
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);
k3 = omega*sqrt(mu0*ep0*ep3);
k4 = omega*sqrt(mu0*ep0*ep4);




% Source Location
zp = 5e-3;

% Layer Heights
z1 = 0;
z2 = 4e-3;
z3 = 6e-3;
z4 = 10e-3;

% Layer thickness
d1 = z2 - z1;
d2 = z3 - z2;
d3 = z4 - z3;
dp = zp - z2;
dpp = z3 - zp;
% TE/TM switch
nu = 1;

kz1 =  sqrt(k1 ^2 - kp .^2);

kz2 =  sqrt(k2 ^2 - kp .^2);

kz3 =  sqrt(k3 ^2 - kp .^2);

kz4 =  sqrt(k4 ^2 - kp .^2);

% Define impedances
if nu == 0
    
    % TE Case
    Z1 =  omega./kz1;
    Z2 =  omega./kz2;
    Z3 =  omega./kz3;
    Z4 =  omega./kz4;
else
    % TM case
    Z1 =  kz1./(omega*ep1);
    Z2 =  kz2./(omega*ep2);
    Z3 =  kz3./(omega*ep3);
    Z4 =  kz4./(omega*ep4);
end


Gamma_12 =  (Z1 - Z2)./(Z1 + Z2);
Gamma_32 =  (Z3 - Z2)./(Z3 + Z2);
Gamma_43 =  (Z4 - Z3)./(Z4 + Z3);


% % Reflection Coefficients
Gamma_left =  (Gamma_12 + (-1)*exp(-2i*kz1*d1))...
    ./(1 + Gamma_12*(-1).*exp(-2i*kz1*d1)); % Left-looking
Gamma_right =  (Gamma_32 + Gamma_43.*exp(-2i*kz3*d3))...
    ./(1 + Gamma_32.*Gamma_43.*exp(-2i*kz3*d3)); % Right-looking

% Unknown A
A =  (Gamma_left .* exp(1i*kz2*2*z2))./(1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * d2)) ...
    .*( exp(-1i * kz2 * zp) + Gamma_right .* exp( -1i * kz2 * (2*z3 - zp)));

% Unknown B
B =  (Gamma_right .* exp(-1i*kz2*2*z3))./(1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * d2)) ...
    .*( exp(+1i * kz2 * zp) + Gamma_left .* exp( +1i * kz2 * (2*z2 - zp)));

% Denominator
D =  1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * dpp).*exp(-2i * kz2 * dp);


end
