function D = TLGFd(kp)
% ***********************************************************************
%
%      Computes the Transmission Line Greens Function of a three
%      layer structure with source located at the center of the middle
%      layer
%
% ***********************************************************************

% Data from datafile of test_case3 from the reference paper
f = 1e12;
omega = 2*pi*f;
lambda = 3e8/f;

% Example Validations

% Material Properties
ep1 = 1; % Air
ep2 = 9.7 ; % GaN/AlGaN layers combined
ep3 = 1; % Silicon base


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

kz1 =  sqrt(k1 ^2 - kp .^2);

kz2 =  sqrt(k2 ^2 - kp .^2);

kz3 =  sqrt(k3 ^2 - kp .^2);

% Enforce kzn on the top sheet
for i = 1 : length(kz1)
    if imag(kz1(i)) <= 0
        %         kz1(i) = conj(kz1(i));
        kz1(i) = -(kz1(i));
    end
    
%     if imag(kz2(i)) <= 0
%         kz2(i) = conj(kz2(i));
%     end
    
    if imag(kz3(i)) <= 0
        %         kz3(i) = conj(kz3(i));
        kz3(i) = -(kz3(i));
    end
end

% Define impedances
if nu == 0
    
    % TE Case
    Z1 =  omega./kz1;
    Z2 =  omega./kz2;
    Z3 =  omega./kz3;
else
    % TM case
    Z1 =  kz1./(omega*ep1);
    Z2 =  kz2./(omega*ep2);
    Z3 =  kz3./(omega*ep3);
end

% Normalize Impedances to free-space
% Z0 = sqrt(mu0/ep0);
% Z1 = Z1/Z0;
% Z2 = Z2/Z0;
% Z3 = Z3/Z0;

%% Reflection Coefficients
Gamma_left = (Z3 - Z2) ./ (Z3 + Z2); % Left-looking
Gamma_right =  (Z1 - Z2) ./ (Z1 + Z2); % Right-looking

% Unknown A
% A =  (Gamma_left .* exp(1i*kz2*2*z2))./(1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * d2)) ...
%     .*( exp(-1i * kz2 * zp) + Gamma_right .* exp( -1i * kz2 * (2*z3 - zp)));

% Unknown B
% B =  (Gamma_right .* exp(-1i*kz2*2*z3))./(1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * d2)) ...
%     .*( exp(+1i * kz2 * zp) + Gamma_left .* exp( +1i * kz2 * (2*z2 - zp)));

% Denominator
D =  1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * d);

end